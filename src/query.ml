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

(* query a fragment database using
   - a query fragment
   - a maximum allowed RMSD to this query *)

open Core
open Printf
open Indexed_frag
open Frag_index
open Log

module DB      = Frag_db
module FragSet = Frag_index.FragSet
module HT      = Caml.Hashtbl
module Idx     = Frag_index
module L       = List
module MU      = My_utils
module S       = String

type option_type =
    { index     : string ref ;
      idx       : string ref ;
      db        : string ref ;
      query     : string ref ;
      out_f     : string ref ;
      filt      : string ref ;
      mode      : string ref ;
      io_mode   : DB.io_mode ref;
      d         : float  ref ;
      dx        : float  ref ;
      rand_by   : (int * int) ref;
      min       : int    ref ;
      max       : int    ref ;
      n         : int    ref ;
      excl      : string ref ;
      res       : string ref ;
      np        : int    ref ;
      no_coords : bool   ref ;
      force     : bool   ref ;
      debug     : bool   ref }

let random_residue res_num =
  let rand_f () =
    (* four digits before decimal point is the max allowed by PDB format
       for coordinates *)
    Random.float 99.
  in
  sprintf
  "ATOM  %5d  N   ALA A%4d    %8.3f%8.3f%8.3f  1.00 30.00      A    N  \n\
   ATOM  %5d  CA  ALA A%4d    %8.3f%8.3f%8.3f  1.00 30.00      A    C  \n\
   ATOM  %5d  C   ALA A%4d    %8.3f%8.3f%8.3f  1.00 30.00      A    C  \n\
   ATOM  %5d  O   ALA A%4d    %8.3f%8.3f%8.3f  1.00 30.00      A    O  \n"
    (4 * (res_num - 1) + 1) res_num (rand_f()) (rand_f()) (rand_f())
    (4 * (res_num - 1) + 2) res_num (rand_f()) (rand_f()) (rand_f())
    (4 * (res_num - 1) + 3) res_num (rand_f()) (rand_f()) (rand_f())
    (4 * (res_num - 1) + 4) res_num (rand_f()) (rand_f()) (rand_f())
(* toplevel test:
   printf "%s" (random_residue 1);;
*)

type user_choice = Res_list of string
                 | Nb_res of int

let res_list_or_nb_res res_list nb_res =
  if res_list <> "" then
    Res_list res_list
  else
    Nb_res nb_res

let string_of_user_choice = function
  | Res_list rl -> rl
  | Nb_res n -> string_of_int n

let fuzzy_query ranker res query d force index =
  info "index: %s" index;
  (* distance between reference and query *)
  let ref_frag_fn    = Filename.temp_file "ref_" ".pdb" in
  let ref_frag_lines = MU.read_while ((<>) "END") index in
  let nb_res = ((L.length ref_frag_lines) - 1) / 4 in
  MU.string_list_to_file ref_frag_fn ref_frag_lines;
  let d_rq =
    fst
      (Idx.rank_id_of_string
         (MU.get_command_output
            (sprintf "%s %s %s %s 0"
               ranker
               (string_of_user_choice (res_list_or_nb_res res nb_res))
               ref_frag_fn
               query))) in
  Sys.remove ref_frag_fn;
  let d_inf = d_rq -. d in
  let d_sup = d_rq +. d in
  (* load fragments index to ref *)
  let some_fragments = Idx.load_index index force in
  (* output only the interesting subset of the fragments DB *)
  let _, ge_d_inf = FragSet.split_ge (d_inf, []) some_fragments in
  let potential_frags, _ = FragSet.split_le (d_sup, []) ge_d_inf in
  (* this set will be reduced later using results from queries using
     different references *)
  info "first: %f" (fst (FragSet.min_elt potential_frags));
  info "last: %f" (fst (FragSet.max_elt potential_frags));
  let requested = HT.create 10_000_000 in
  FragSet.iter
    (fun (_rmsd, frag_ids) ->
       L.iter frag_ids
         ~f:(fun frag_id -> HT.add requested frag_id ()))
    potential_frags;
  info "nb candidates: %d" (HT.length requested);
  (nb_res, requested)

let read_fst_frag fn =
  let frags_seen = ref 0 in
  let still_fst_frag l =
    if MU.starts_with l "TER "
    then (MU.plus_plus frags_seen) = 1
    else true
  in
  MU.read_while still_fst_frag fn

let create_index opts ranker ref_frag all_fragments_fn out_fn =
  let ref_frag_fn = Filename.temp_file "rand_ref_" ".pdb" in
  let db_index = DB.load_db_index !(opts.db) false in
  let frag_db_fd = Unix.openfile [Unix.O_RDONLY] !(opts.db) in
  let fst_frag_lines = read_fst_frag all_fragments_fn in
  let nb_res = ((L.length fst_frag_lines) - 1) / 4 in
  MU.with_out_file
    ref_frag_fn
    (DB.output_several_frags
       false DB.Txt2txt no_sequence_filter db_index frag_db_fd nb_res [ref_frag]);
  Unix.close frag_db_fd;
  info "creating %d residues RMSD index for %s in %s..."
    nb_res ref_frag ref_frag_fn;
  Random.self_init(); (* varies each run *)
  (* dump ref_frag as the header in the index file too *)
  MU.run_command (
    sprintf "cat %s > %s; echo END >> %s" ref_frag_fn out_fn out_fn);
  (* append all RMSDs in there too *)
  MU.run_command
    (sprintf "%s %d %s %s 0 >> %s"
       ranker
       nb_res
       all_fragments_fn
       ref_frag_fn
       out_fn
    );
  Sys.remove ref_frag_fn

let merge fquery_results =
  L.reduce_exn fquery_results
    ~f:(fun (nb_res1, frags1) (nb_res2, frags2) ->
          if nb_res1 = nb_res2 then
            let new_set = MU.inter frags1 frags2 in
            let _ = info "nb frags: %d" (HT.length new_set) in
            (nb_res1, new_set)
          else failwith
            (sprintf
               ("indexes have different fragment sizes: %d Vs. %d")
               nb_res1 nb_res2))

let exact_query opts ranker db_index frag_db_fd nb_res requested =
  (* EXACT QUERY *)
  info "loading %d fragments ..." (L.length requested);
  let seq_filter = sequence_filter !(opts.filt) in
  MU.with_out_file !(opts.out_f)
    (DB.output_several_frags
       false !(opts.io_mode) seq_filter db_index frag_db_fd nb_res requested);
  (* only keep the ones within query distance to the fragment query *)
  info "filtering them ...";
  let choice = res_list_or_nb_res !(opts.res) nb_res in
  let optim = match choice with
    | Res_list _ -> ""
    | Nb_res _ -> " -bin" in
  let rmsds_str =
    MU.get_command_output
      (sprintf "%s %s %s %s %f%s"
         ranker
         (string_of_user_choice choice)
         !(opts.out_f)
         !(opts.query)
         !(opts.d)
         optim) in
  info "writing them ...";
  let nb_survivors = ref 0 in
  let survived =
    L.fold (S.split rmsds_str ~on:'\n') ~init:[]
      ~f:(fun acc rmsd_frag ->
            if rmsd_frag = "" then acc
            else match S.split rmsd_frag ~on:' ' with
              | [rmsd; id] ->
                let rmsd_f = Float.of_string rmsd in
                incr nb_survivors;
                (rmsd_f, id) :: acc
              | _ ->
                let _ = error "unexpected line in the output of ranker: %s"
                    rmsd_frag in
                acc) in
  info "Found %d fragment(s)" !nb_survivors;
  (nb_res, survived, !nb_survivors)

let merge_query_outputs
    (nb_res1, survived1, nb_survivors1)
    (nb_res2, survived2, nb_survivors2) =
  assert (nb_res1 = nb_res2);
  (nb_res1,
   L.append survived1 survived2,
   nb_survivors1 + nb_survivors2)

let fuzzy_then_exact_query opts ranker db_index frag_db_fd =
  (* FUZZY QUERY *)
  let indexes = S.split ~on:',' !(opts.idx) in
  let fquery_results =
    if !(opts.np) > 1 then
      Parmap.parmap ~ncores:!(opts.np)
        (fuzzy_query
           ranker !(opts.res) !(opts.query) !(opts.d) !(opts.force))
        (Parmap.L indexes)
    else
      L.map indexes
        ~f:(fuzzy_query
              ranker !(opts.res) !(opts.query) !(opts.d) !(opts.force))
  in
  let nb_res, requested_set = merge fquery_results in
  let requested =
    let to_filter_out = MU.string_set_of_file !(opts.excl) in
    let lst = MU.keys requested_set in
    L.filter lst
      ~f:(fun id ->
            let pdb, _chain, _maybe_id = deconstruct_id id in
            not (String.Set.mem to_filter_out pdb))
  in
  let n, by = !(opts.rand_by) in
  if n = 0 then
    exact_query opts ranker db_index frag_db_fd nb_res requested
  else
    let to_process = MU.split_by by requested in
    let rec loop ((_, _, nb_found) as acc) = function
      | [] -> acc
      | x :: xs ->
        if nb_found >= n then
          acc
        else
          let partial = exact_query opts ranker db_index frag_db_fd nb_res x in
          loop (merge_query_outputs acc partial) xs
    in loop (nb_res, [], 0) to_process

let parse_min_option_string opts s =
  if s = "" then ()
  else
    let n, dx =
      try Scanf.sscanf s "%d:%f" (fun n dx -> (n, dx))
      with _ -> failwith ("query.ml: parse_min_option_string: \
                           cannot parse: " ^ s)
    in
    opts.min := n;
    opts.dx := dx

let parse_max_option_string opts s =
  if s = "" then ()
  else
    let n, dx =
      try Scanf.sscanf s "%d:%f" (fun n dx -> (n, dx))
      with _ -> failwith ("query.ml: parse_max_option_string: \
                           cannot parse: " ^ s)
    in
    opts.max := n;
    opts.dx := -.dx

let parse_rand_by_option_string opts s =
  if s = "" then ()
  else
    opts.rand_by :=
      try Scanf.sscanf s "%d:%d" (fun i j -> (i, j))
      with _ -> failwith ("query.ml: parse_rand_by_option_string: \
                           cannot parse: " ^ s)

let parse_index_option_string s =
  match S.split s ~on:':' with
    | [ref_frag; all_frags] -> (ref_frag, all_frags)
    | _ -> failwith ("query.ml: parse_index_option_string: cannot parse: " ^ s)

let main () =

  set_log_level INFO;
  color_on();

  let opts = {
    index     = ref ""          ;
    idx       = ref ""          ;
    db        = ref ""          ;
    query     = ref ""          ;
    out_f     = ref ""          ;
    filt      = ref ""          ;
    mode      = ref "t2t"       ;
    io_mode   = ref DB.Txt2txt  ;
    d         = ref 0.0         ;
    dx        = ref 0.1         ;
    rand_by   = ref (0, 0)      ;
    min       = ref 0           ;
    max       = ref 0           ;
    n         = ref 0           ;
    res       = ref ""          ;
    np        = ref 1           ;
    no_coords = ref false       ;
    excl      = ref "/dev/null" ;
    force     = ref false       ;
    debug     = ref false       } in
  let cmd_line = Arg.align [
    "-d"  , Arg.Set_float  opts.d    , "distance_to_query (> 0.0)"       ;
    "-db" , Arg.Set_string opts.db   , "fragments_DB_file"               ;
    "-f"  , Arg.Set        opts.force, " force reindex (slow)"           ;
    "-idx", Arg.Set_string opts.idx  , "a,b,c coma-separated list of \
                                        RMSD index files"                ;
    "-i"  , Arg.Set_string opts.index, "ref_frag_id:all_fragments_file \
                                        create a new RMSD index"         ;
    "-n"  , Arg.Set_int    opts.n    , "n_best only output the n \
                                        best fragments"                  ;
    "-m"  , Arg.Set_string opts.mode , " [MODE] mode may be one of \
        {t2t|b2t|b2b} \
        (t2t: read and write pdb ATOM lines (default); \
         b2t: read binary coordinates and output pdb ATOM lines; \
         b2b: read and write binary coordinates)"                        ;
    "-x"  , Arg.Set_string opts.excl , "exclude_file one PDB id per line";
    "-min", Arg.String (parse_min_option_string opts),
            "N:dx min_nb_frags:query_dist_delta (incompatible with -max)";
    "-max", Arg.String (parse_max_option_string opts),
            "N:dx max_nb_frags:query_dist_delta (incompatible with -min)";
    "-np" , Arg.Set_int    opts.np   , "ncores"                          ;
    "-nc" , Arg.Set        opts.no_coords, " no coordinates written out, \
                                             only sequences"             ;
    "-o"  , Arg.Set_string opts.out_f, "output_file"                     ;
    "-q"  , Arg.Set_string opts.query, "query_fragment"                  ;
    "-res", Arg.Set_string opts.res  , "1,2,8,9 coma-separated list \
                                        of residues"                     ;
    "-rby", Arg.String (parse_rand_by_option_string opts),
            "r:N min_nb_frags:group_size"                                ;
    "-seq", Arg.Set_string opts.filt , " sequence_regexp \
                                         (use 1 char AA codes)"          ;
    "-v"  , Arg.Set        opts.debug, " debug mode"                     ] in
  let use_msg = sprintf
    "usage: %s -db fragments_file -q query_frag -d 0.5 -n 10 -o out \
               -idx fragments_RMSD_index [-x exclude_file] [-f] [-v] \
                                         [-m {t2t|b2t|b2b}]"
    Sys.argv.(0) in
  Arg.parse cmd_line ignore use_msg;
  let indexing = !(opts.index) <> "" in
  if (not indexing && !(opts.idx)   = "")  ||  (* no RMSD index *)
     (not indexing && !(opts.db)    = "")  ||  (* no DB *)
     (not indexing && !(opts.query) = "")  ||  (* no query *)
     (not indexing && !(opts.d)    <= 0.0) ||  (* no RMSD tolerance *)
                      !(opts.out_f) = ""   ||  (* no output file *)
    ((!(opts.min) <> 0) && (!(opts.max) <> 0)) (* incompatible options *)
  then begin
    fatal "%s" use_msg;
    exit 1
  end;
  let ranker =
    if MU.detect_cmd "ranker_aa" then "ranker_aa"
    else MU.getenv_or_fail "RANKER"
      "the RANKER environment variable must point \
       to ranker_aa exe or the ranker_aa command must be \
       in your PATH" in
  if indexing then
    (let ref_frag, all_frags = parse_index_option_string !(opts.index) in
     create_index opts ranker ref_frag all_frags !(opts.out_f);
     exit 0);
  opts.io_mode := DB.io_mode_of_string !(opts.mode);
  let db_index = DB.load_db_index !(opts.db) !(opts.debug) in
  let frag_db_fd = Unix.openfile [Unix.O_RDONLY] !(opts.db) in
  let query_output =
    ref (fuzzy_then_exact_query opts ranker db_index frag_db_fd) in
  (* -min or -max option *)
  while ((!(opts.min) <> 0) && ((trd3 !query_output) < !(opts.min))) ||
        ((!(opts.max) <> 0) && ((trd3 !query_output) > !(opts.max))) do
    (* update query distance *)
    opts.d := !(opts.d) +. !(opts.dx);
    info "query distance updated: %f" !(opts.d);
    query_output := fuzzy_then_exact_query opts ranker db_index frag_db_fd;
  done;
  let nb_res, survived, _nb_found = !query_output in
  let selected =
    L.sort survived
      ~cmp:(fun (rmsd1, _) (rmsd2, _) -> Float.ascending rmsd1 rmsd2) in
  MU.with_out_file !(opts.out_f)
    (fun out ->
       let n = !(opts.n) in
       let lst =
         if n <> 0 then
           let _ = info "Will output %d" n in
           L.take selected n
         else
           selected
       in
       L.iter lst
         ~f:(fun (rmsd, id) ->
               fprintf out "REMARK %f\n" rmsd;
               (* FBR: less efficent than requesting many at a time *)
               DB.output_one_frag !(opts.no_coords) !(opts.io_mode)
                 no_sequence_filter db_index frag_db_fd nb_res id out));
  info "Done";
  Unix.close frag_db_fd
;;

main()
