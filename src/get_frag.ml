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

(* get some fragments out of a fragments DB, either randomly
   or selected by id. *)

open Core.Std
open Printf
open Indexed_frag
open Log

module DB = Frag_db
module HT = Caml.Hashtbl
module L  = List
module MU = My_utils

let int_of_bool = function
  | true  -> 1
  | false -> 0

let main () =

  set_log_level INFO;
  color_on();

  let db        = ref ""    in
  let in_f      = ref ""    in
  let frag_id   = ref ""    in
  let out_f     = ref ""    in
  let filt      = ref ""    in
  let mode      = ref "t2t" in
  let n         = ref 0     in
  let no_coords = ref false in
  let r         = ref 0     in
  let debug     = ref false in
  let cmd_line = Arg.align
    [ "-db" , Arg.Set_string db     , " pdb_chains_file"                      ;
      "-i"  , Arg.Set_string in_f   , " input_file (one fragment id per line, \
                                       incompatible with -r)"                 ;
      "-id" , Arg.Set_string frag_id, " fragment_id (incompatible w/ -i/-r)"  ;
      "-m"  , Arg.Set_string mode   , " [MODE] mode may be one of \
        {t2t|b2t|b2b} \
        (t2t: read and write pdb ATOM lines (default); \
         b2t: read binary coordinates and output pdb ATOM lines; \
         b2b: read and write binary coordinates)"                             ;
      "-n"  , Arg.Set_int    n      , " fragment_length (> 0)"                ;
      "-nc" , Arg.Set no_coords     , " no coordinates written out, \
                                        only sequences"                       ;
      "-o"  , Arg.Set_string out_f  , " output_file (default is stdout)"      ;
      "-r"  , Arg.Set_int    r      , " nb_random_fragments \
                                       (incompatible with -i)"                ;
      "-seq", Arg.Set_string filt   , " sequence_regexp (use one char \
                                        amino acid codes)"                    ;
      "-v"  , Arg.Set        debug  , " debug mode"                           ]
  in
  let use_msg = sprintf
    "usage: %s -db pdb_chains_file\n\
              {-id frag_id | -i fragments_id_file | -r N -n L}\n\
              [-o out] [-v] [-help]"
    Sys.argv.(0) in
  Arg.parse cmd_line ignore use_msg;
  if !db = "" || (* no DB file given *)
     (((int_of_bool (!in_f <> "")) +
       (int_of_bool (!r > 0)) +
       (int_of_bool (!frag_id <> ""))) > 1) (* incompatible options *)
  then begin
    fatal (lazy use_msg);
    exit 1
  end;
  let index = DB.load_db_index !db !debug in
  let frag_db_fd = Unix.openfile [Unix.O_RDONLY] !db in
  let requested =
    if !frag_id <> "" then [!frag_id]
    else if !in_f <> "" then In_channel.read_lines !in_f
    else if !r > 0 then
      let long_enough_sequences =
        MU.keys index
        (* L.filter *)
        (*   (MU.keys index) *)
        (*   (fun key -> *)
        (*     let value = HT.find index key in *)
        (*     let len = String.length value.sequence in *)
        (*     len <= !n) *)
      in
      let rand_pdbs =
        L.permute
          (* varies each time the software is run *)
          ~random_state:(Random.State.make_self_init())
          long_enough_sequences
      in
      let r_rand_pdbs = L.take rand_pdbs !r in
      (* fragment index will be chosen randomly *)
      L.map r_rand_pdbs (fun s -> s ^ "_~")
    else []
  in
  let seq_filter = sequence_filter !filt in
  let io_mode = DB.io_mode_of_string !mode in
  DB.reset_count();
  MU.with_out_file
    (if !out_f <> "" then !out_f else "/dev/stdout")
    (DB.output_several_frags
       !no_coords io_mode seq_filter index frag_db_fd !n requested);
  Unix.close frag_db_fd;
  let count = DB.get_count() in
  info (lazy (sprintf "wrote %d fragments of %d residues" count !n))
;;

main()
