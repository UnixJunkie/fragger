(* Copyright (C) 2013, Zhang Initiative Research Unit,
 * Advance Science Institute, Riken
 * 2-1 Hirosawa, Wako, Saitama 351-0198, Japan
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License, with
 * the special exception on linking described in file LICENSE.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. *)

(* Keep only quads of lines N,CA,C,O per residue from the input PDB *)

open Core.Std
open Indexed_frag
open Printf
open Log

module AA   = Amino_acid
module Buff = Buffer
module F    = Filename
module HT   = Caml.Hashtbl
module L    = List
module MU   = My_utils
module PDB  = Pdb_parser
module S    = String

(* return (only_first_chain_lines, remaining_lines)
   n.b. aborts early, not the standard partition function for lists *)
let partition_by_chain pdb_lines =
  match pdb_lines with
      []     -> ([], [])
    | l :: _ ->
        let id = PDB.chain_id_of_pdb_line l in
        let rec loop lst acc = match lst with
            []       -> (L.rev acc, [])
          | l' :: ls ->
              if Char.equal id (PDB.chain_id_of_pdb_line l') then
                loop ls (l' :: acc)
              else
                (L.rev acc, lst)
        in
        loop pdb_lines []

(* return (only_1st_residue_lines, remaining_lines)
   n.b. aborts early, not the standard partition function for lists *)
let partition_by_residue pdb_lines =
  match pdb_lines with
      []     -> ([], [])
    | l :: _ ->
        let id = PDB.res_num_of_pdb_line l in
        let rec loop lst acc = match lst with
            []       -> (L.rev acc, [])
          | l' :: ls ->
              if Int.equal id (PDB.res_num_of_pdb_line l') then
                loop ls (l' :: acc)
              else
                (L.rev acc, lst)
        in
        loop pdb_lines []

(* return [residue_1_lines; residue_2_lines; ...] *)
let group_by_residue pdb_lines =
  MU.group_by partition_by_residue pdb_lines

(* detect holes in residue numbering *)
let rec detect_chain_break = function
  | res1 :: res2 :: others ->
    (try
       let n1 = PDB.res_num_of_pdb_line (L.hd_exn res1) in
       let n2 = PDB.res_num_of_pdb_line (L.hd_exn res2) in
       if n2 <> n1 + 1 then true (* chain break *)
       else (* examine further down *)
         detect_chain_break (res2 :: others)
     with _ -> true) (* impossible to read resnum considered a chain break *)
  | _ -> false

(* return [chain_A_lines; chain_B_lines; ...] *)
let group_by_chain pdb_lines =
  MU.group_by partition_by_chain pdb_lines

let only_NCACO_quads pdb_name warn1 warn2 l =
  L.filter l
    ~f:(fun residue ->
          match residue with
              [n; ca; c; o] ->
                if PDB.bb_atom_name_of_pdb_line n  = "N " &&
                   PDB.bb_atom_name_of_pdb_line ca = "CA" &&
                   PDB.bb_atom_name_of_pdb_line c  = "C " &&
                   PDB.bb_atom_name_of_pdb_line o  = "O "
                then
                  true
                else
                  let _ =
                    if !warn1 then
                      (warn
                         (lazy (sprintf "bad bb atom types in %s" pdb_name));
                       warn1 := false);
                  in
                  false
            | _ ->
              let _ =
                if !warn2 then
                  (warn
                     (lazy (sprintf "bad number of bb atoms in %s" pdb_name));
                   warn2 := false);
              in
              false)

(* string of one letter amino acids for all quads of a chain *)
let sequence chain_atom_lines =
  let buff = Buff.create 1000 in
  L.iter chain_atom_lines
    ~f:(fun residue ->
          let aa_str = PDB.res_name_of_pdb_line (L.hd_exn residue) in
          let aa = AA.aa_one_of_string aa_str in
          Buff.add_char buff (AA.char_of_aa_one aa));
  Buffer.contents buff

(* write a pdb_chain to the database and return its descriptor *)
let write_pdb_chain offset out pdb chain seq quads =
  let pdb_chain = sprintf "%s_%c" pdb chain in
  let header = sprintf "TER %s %s\n" pdb_chain seq in
  let ter_line_len = S.length header in
  (* TER line *)
  fprintf out "%s" header;
  (* ATOM lines *)
  let atom_lines_len = ref 0 in
  L.iter quads
    ~f:(L.iter
          ~f:(fun l ->
                (*                                + 1 for '\n' *)
                atom_lines_len := !atom_lines_len + 1 + (S.length l);
                fprintf out "%s\n" l));
  (* ATOM xyz coordinates in binary format *)
  let buff = Buffer.create 8000 in
  L.iter quads
    ~f:(L.iter
          ~f:(fun atom_line ->
                let x, y, z = PDB.xyz_of_pdb_line atom_line in
                Buffer.add_string buff (MU.str_of_float32 x);
                Buffer.add_string buff (MU.str_of_float32 y);
                Buffer.add_string buff (MU.str_of_float32 z)));
  { pdb_chain      = pdb_chain            ;
    offset         = !offset              ;
    ter_line_len   = ter_line_len         ;
    atom_lines_len = !atom_lines_len      ;
    sequence       = seq                  ;
    bin_coords     = Buffer.contents buff }

let process_chain no_dup_seq seen_sequences index out offset pdb_in
    warn1 warn2 pdb _secondary lines =
  match lines with
  | [] -> ()
  | l :: _ ->
    let chain = PDB.chain_id_of_pdb_line l in
    try
      let residues = group_by_residue lines in
      if detect_chain_break residues
      then raise (Failure "broken chain");
      let quads = only_NCACO_quads pdb_in warn1 warn2 residues in
      let seq = sequence quads in
      if no_dup_seq && (HT.mem seen_sequences seq) then ()
      else begin
        HT.add seen_sequences seq ();
        let frag_desc = write_pdb_chain offset out pdb chain seq quads in
        let frag_len = block_length frag_desc in
        offset := Int64.(!offset + Int64.of_int frag_len);
        HT.add index frag_desc.pdb_chain frag_desc
      end
    with
      | Failure msg ->
        error (lazy (sprintf "problem while parsing at residue level in %s %c"
                       pdb_in chain));
        error (lazy msg)
      | AA.Unknown_aa_one aa ->
        warn (lazy (sprintf "non standard AA in %s %c: %s" pdb_in chain aa))

let fragment_pdb no_dup_seq index out offset pdb_in =
  let pdb = F.chop_extension (F.basename pdb_in) in
  let lines = MU.some_lines_of_file (PDB.is_atom_no_altloc) pdb_in in
  let secondary = [] in
  (* PDB.parse_ksdssp_output (PDB.run_ksdssp ksdssp pdb_in) in *)
  let chains =
    try group_by_chain lines
    with Failure msg ->
      error (lazy (sprintf "problem while parsing at chain level in %s"
                     pdb_in));
      error (lazy msg);
      [] (* PDB is ignored if error at this level of parsing *)
  in
  let warn1 = ref true in
  let warn2 = ref true in
  let seen_sequences = HT.create 30 in
  L.iter chains
    ~f:(process_chain
          no_dup_seq
          seen_sequences
          index
          out
          offset
          pdb_in
          warn1
          warn2
          pdb
          secondary)

let main () =

  set_log_level INFO;
  color_on();

  let pdb_in     = ref ""    in
  let out_f      = ref ""    in
  let np         = ref 1     in
  let no_dup_seq = ref false in
  let cmd_line =
    Arg.align ["-i"     , Arg.Set_string pdb_in, "bb_pdb_files_list";
               "-np"    , Arg.Set_int np       , "ncores"           ;
               "-nodups", Arg.Set no_dup_seq   , " rm duplicated \
                          sequences from a PDB"                     ;
               "-o"     , Arg.Set_string out_f , "out_file"         ] in
  let usage_msg =
    sprintf
      "usage: %s -i pdb_files_list -o out_file [-np ncores]"
      Sys.argv.(0) in
  Arg.parse cmd_line ignore usage_msg;
  if !out_f = "" || !pdb_in = "" then begin
    fatal (lazy usage_msg);
    exit 1;
  end;
  let pdbs = MU.string_list_of_file !pdb_in in
  let process_pdb = (fun pdb -> (pdb, [])) in
  let pdb_ss =
    if !np > 1 then
      Parmap.parmap
        ~ncores:!np
        process_pdb
        (Parmap.L pdbs)
    else
      List.map pdbs ~f:process_pdb
  in
  let index = HT.create 50_000 in
  let output_file = !out_f in
  let idx_file = !out_f ^ ".idx" in
  MU.with_out_file output_file
    (fun out ->
      let offset = ref 0L in
      (L.iter pdb_ss
         ~f:(fun (pdb, _ss) ->
               fragment_pdb !no_dup_seq index out offset pdb)));
  info (lazy "saving DB index for next time ...");
  MU.with_out_file ~bin:true idx_file
    (fun out -> Marshal.to_channel out index [Marshal.No_sharing]);
  info (lazy "saved")
;;

main()
