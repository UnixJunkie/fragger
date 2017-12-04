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

(* to create/load the RMSD index of many fragments
   to a given reference fragment *)

open Core
open Printf
open Log

module A  = Array
module L  = List
module MU = My_utils

module FragSet = struct
  include Set_with_find.Make(struct
    type t = (float * string list)
    let compare (rmsd1, _ids1) (rmsd2, _ids2) = Float.compare rmsd1 rmsd2
  end)
  let of_list l = L.fold l ~init:empty ~f:(fun acc x -> add x acc)
  (* return two sets: (elts < elt, elts >= elt) *)
  let split_ge elt s =
    let l, present, r = split elt s in
    if present then
      l, add (find elt s) r
    else
      l, r
  (* return two sets: (elts <= elt, elts > elt) *)
  let split_le elt s =
    let l, present, r = split elt s in
    if present then
      add (find elt s) l, r
    else
      l, r
end

let rank_id_of_string s =
  try Scanf.sscanf s "%f %s" (fun dist id -> (dist, id))
  with _ ->
    failwith ("frag_index.ml: rank_id_of_string: unable to parse: " ^ s)

(* group fragments with same distance to ref
   PRECONDITION: lst must be sorted by dist to ref *)
let group_by_dist lst =
  let rec group_by_dist_priv l acc =
    match l with
      | [] -> acc
      | (dist1, id1) :: xs ->
        match acc with
          | [] -> group_by_dist_priv xs [(dist1, [id1])]
          | (dist2, ids) :: ys ->
            if Float.compare dist1 dist2 = 0 then
              group_by_dist_priv xs ((dist1, id1 :: ids) :: ys)
            else
              group_by_dist_priv xs ((dist1, [id1]) :: acc)
  in group_by_dist_priv lst []

let cmp_rank_id (rmsd1, _) (rmsd2, _) =
  Float.compare rmsd1 rmsd2

let string_of_rank_id (rank, id) =
  sprintf "%f %s" rank id

let load_index idx_file force =
  let bin_idx_file = idx_file ^ ".bin" in
  if (Sys.file_exists bin_idx_file = `Yes) && (not force) then begin
    info "loading index from previous run ...";
    let fragments =
      In_channel.with_file ~binary:true bin_idx_file
        ~f:(fun input -> Marshal.from_channel input) in
    info "loaded";
    fragments
  end else begin
    info "loading index from text file ...";
    let before = Unix.gettimeofday() in
    let nb_frags = ref 0 in
    let parsed_lines =
      In_channel.with_file idx_file
        ~f:(fun input ->
              (* skip the reference fragment *)
              while Caml.input_line input <> "END" do () done;
              In_channel.fold_lines input
                ~init:[]
                ~f:(fun acc l ->
                      let _ = incr nb_frags in
                      (rank_id_of_string l) :: acc)) in
    let sorted_fragments = L.sort parsed_lines ~cmp:cmp_rank_id in
    let grouped_fragments = group_by_dist sorted_fragments in
    let fragments = FragSet.of_list grouped_fragments in
    let after = Unix.gettimeofday() in
    info "loaded index for %d fragments in %.3fs"
      !nb_frags (after -. before);
    info "saving index for next time ...";
    MU.with_out_file ~bin:true bin_idx_file
      (fun out -> Marshal.to_channel out fragments [Marshal.No_sharing]);
    info "saved";
    fragments
  end
