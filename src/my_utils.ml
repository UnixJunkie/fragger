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

(* to work in the toplevel
   #use "topfind";;
   #thread;;
   #require "core.top";;
   #require "str";;
   #use "my_utils.ml";;
*)

open Core

module A   = Array
module F   = Filename
module HT  = Caml.Hashtbl
module L   = List
module BL  = BatList
module Mar = Marshal
module P   = Printf
module RNG = Random
module S   = String

(* ============================== Operators ===================== *)

(* function composition operator (chaining functions) *)
let ($@) f g x = f (g x)
(* in case the previous operator give problems while parsing *)
let compose f g x = f (g x)

let (|>) x f = f x

let ignore_fst _ y = y
let ignore_snd x _ = x

let tap f x =
  f x;
  x

let guarded_tap p f x =
  if p then f x;
  x

let n_times n f =
  for i = 1 to n do
    f()
  done

let plus_plus ri =
  incr ri;
  !ri

(* ============================== Conversion functions ===================== *)

(* ASCII to integer *)
let atoi = Int.of_string

(* integer to ASCII *)
let itoa = string_of_int

let atof = Float.of_string

let ftoa = string_of_float

let itof = Float.of_int

(* truncate, does not round *)
let ftoi = Int.of_float

(* ============================== I/O ============================== *)

(* like Python's readlines() *)
let string_list_of_file f =
  In_channel.read_lines f

(* store a list into a file optionally adding header/footer *)
let list_to_file ?(head = "") ?(foot = "") f to_string l =
  let out = open_out f in
  if head <> "" then P.fprintf out "%s" head;
  L.iter
    ~f:(fun elt -> P.fprintf out "%s" (to_string elt))
    l;
  if foot <> "" then P.fprintf out "%s" foot;
  close_out_noerr out

(* output a string list to a file, '\n' are added
   optionally adding header/footer *)
let string_list_to_file ?(head = "") ?(foot = "") f l =
  list_to_file ~head:head ~foot:foot f
    (fun s -> P.sprintf "%s\n" s)
    l

(* read fn as long as (p current_line) is true *)
let read_while p fn =
  let res = ref [] in
  In_channel.with_file ~binary:false fn
    ~f:(fun input ->
          let line = ref (Caml.input_line input) in
          while p !line do
            res  := !line :: !res;
            line := Caml.input_line input
          done);
  L.rev !res

(* keep only lines satisfying p *)
let some_lines_of_file p f =
  L.filter ~f:p (string_list_of_file f)

(* call f on lines of file *)
let iter_on_lines_of_file f file =
  let input = open_in file in
  try
    while true do
      f (input_line input)
    done
  with _ -> close_in_noerr input

(* call f on lines of file that satisfy p *)
let iter_on_some_lines_of_file p f file =
  iter_on_lines_of_file
    (fun l -> if p l then f l else ())
    file

(* apply f on lines of file *)
let map_on_lines_of_file f file =
  let res = ref [] in
  iter_on_lines_of_file
    (fun l -> res := (f l) :: !res)
    file;
  L.rev !res

(* apply f on lines of file that satisfy p *)
let map_on_some_lines_of_file p f file =
  let res = ref [] in
  iter_on_lines_of_file
    (fun l ->
       if p l then
         res := (f l) :: !res
       else
         ())
    file;
  L.rev !res

let with_out_file ?(bin = false) fn f =
  let out = (if bin then open_out_bin else open_out) fn in
  let res = f out in
  close_out_noerr out;
  res

(* with_in_file is In_channel.with_file in Janestreet core *)

(* read until the specified number of characters is read *)
let really_read in_fd buf len =
  let was_read = ref 0 in
  while !was_read <> len do
    was_read := !was_read +
      (Unix.read
         ~restart:false ~pos:!was_read in_fd ~buf:buf ~len:(len - !was_read));
  done;
  !was_read

let ftell fd =
  Unix.lseek fd (Int64.of_int 0) Unix.SEEK_CUR

let string_set_of_file fn =
  String.Set.of_list (In_channel.read_lines fn)

let str_of_int32 x = Caml.(
  let high = Int32.to_int (Int32.shift_right_logical x 16) in
  let c1 = Char.chr (high lsr 8) in
  let c2 = Char.chr (high land 0xff) in
  let low = Int32.to_int x in
  let c3 = Char.chr ((low lsr 8) land 0xff) in
  let c4 = Char.chr (low land 0xff) in
  let res = String.create 4 in
  res.[0] <- c1;
  res.[1] <- c2;
  res.[2] <- c3;
  res.[3] <- c4;
  res)

let str_of_float32 x =
  str_of_int32 (Int32.bits_of_float x)

(* read an int32 binary encoded on 4 bytes big-endian *)
let int32_of_str str =
  let c1 = str.[0] in
  let c2 = str.[1] in
  let c3 = str.[2] in
  let c4 = str.[3] in
  let x1 = Caml.Int32.of_int
    (((Caml.Char.code c1) lsl 8) lor (Caml.Char.code c2)) in
  let x2 = Caml.Int32.of_int
    (((Caml.Char.code c3) lsl 8) lor (Caml.Char.code c4)) in
  Caml.Int32.logor (Caml.Int32.shift_left x1 16) x2

(* read a float binary encoded on 4 bytes *)
let float32_of_str str = Caml.(
  Int32.float_of_bits (int32_of_str str)
)
(* test in toplevel:
   float32_of_str (input_line (open_in "test"));;
*)

(* ============================== Strings ============================== *)

(* split 's' using separator 'sep' *)
let str_split sep s =
  Str.split (Str.regexp_string sep) s

(* index of sub string in the super string *)
let str_find ~super ~sub =
  Str.search_forward (Str.regexp_string sub) super 0

(* return all but the first n characters of str *)
let str_tail str ~n =
  S.drop_prefix str n

(* does str starts with pre ? *)
let starts_with str ~pre =
  S.is_prefix str ~prefix:pre

(* ============================== Arrays ============================= *)

let fold_left ~f ~init a =
  let n    = A.length a in
  let acc  = ref init   in
  for i = 0 to n - 1 do
    acc := f !acc a.(i)
  done;
  !acc

(* float array sum *)
let sum_a a = fold_left ~f:(fun a acc -> a +. acc) ~init:0. a

(* float array average *)
let average_a a = (sum_a a) /. (Float.of_int (A.length a))

(* minimum of a float array *)
let min_a a =
  fold_left
    ~f:Float.min
    ~init:Float.max_value
    a

(* maximum of a float array *)
let max_a a =
  fold_left
    ~f:Float.max
    ~init:Float.min_value
    a

(* Bigarray3 with fold *)
module BA3 = struct (* extend type with more operations *)

  include Bigarray.Array3

  let fold arr ~init ~f =
    let acc = ref init in
    let nbi = dim1 arr in
    let nbj = dim2 arr in
    let nbk = dim3 arr in
    for i = 0 to nbi - 1 do
      for j = 0 to nbj - 1 do
        for k = 0 to nbk - 1 do
          acc := f !acc arr.{i, j, k};
        done;
      done;
    done;
    !acc

  let fold_ijk arr ~init ~f =
    let acc = ref init in
    let nbi = dim1 arr in
    let nbj = dim2 arr in
    let nbk = dim3 arr in
    for i = 0 to nbi - 1 do
      for j = 0 to nbj - 1 do
        for k = 0 to nbk - 1 do
          acc := f i j k !acc arr.{i, j, k};
        done;
      done;
    done;
    !acc

  let iter_ijk arr ~f =
    let nbi = dim1 arr in
    let nbj = dim2 arr in
    let nbk = dim3 arr in
    for i = 0 to nbi - 1 do
      for j = 0 to nbj - 1 do
        for k = 0 to nbk - 1 do
          f i j k arr.{i, j, k};
        done;
      done;
    done

  (* blit does a copy, so let's just call it copy *)
  let copy = blit

end;;

(* ============================== Lists ============================== *)

let i_to_j i j =
  let res = ref [] in
  for k = j downto i do
    res := k :: !res
  done;
  !res

let is_empty l =
  match l with
      [] -> true
    | _  -> false

let head l fail_msg = match l with
    [] -> failwith fail_msg
  | hd :: _ -> hd

(* remove first occurence of elt in lst *)
let remove_first elt lst =
  let rec loop l acc = match l with
    | []      -> lst
    | x :: xs ->
        if x = elt
        then L.append (L.rev acc) xs
        else loop xs (x :: acc)
  in
  loop lst []

(* fold_left f on l while p is true *)
let rec fold_while f p acc l =
  match l with
      []      -> acc
    | x :: xs ->
        if p x then
          fold_while f p (f x :: acc) xs
        else
          acc

(* skip all elements of l until p (List.hd l) is true *)
let skip_until p l =
  let rec skip_until_priv lst = match lst with
      []      -> []
    | x :: xs ->
        if p x
        then lst
        else skip_until_priv xs
  in
  skip_until_priv l

(* reversed tail rec. version of stdlib's List.combine *)
let rev_combine l1 l2 =
  let rec combine_priv l1 l2 acc =
    match l1, l2 with
      []       , []        -> acc
    | []       , _
    | _        , []        ->
        failwith "my_utils.ml: rev_combine: list with different lengths"
    | x1 :: xs1, x2 :: xs2 -> combine_priv xs1 xs2 ((x1, x2) :: acc)
  in
  combine_priv l1 l2 []

(* tail rec. version of stdlib's List.combine *)
let combine l1 l2 =
  L.rev (rev_combine l1 l2)

let combine3 l m n =
  let rec loop l1 l2 l3 acc =
    match l1, l2, l3 with
        [], [], [] -> L.rev acc
      | x :: xs, y :: ys, z :: zs ->
          loop xs ys zs ((x, y, z) :: acc)
      | _ -> failwith "my_utils.ml: combine3: different list lengths"
  in
  loop l m n []

let enumerate l =
  let rec enumerate_priv l i acc =
    match l with
        []            -> L.rev acc
      | elt :: others -> enumerate_priv others (i+1) ((i, elt) :: acc)
  in
  enumerate_priv l 0 []

(* prepend (rev 'pre') to 'post' *)
let rec rev_prepend pre post =
  match pre with
    []      -> post
  | x :: xs -> rev_prepend xs (x :: post)

(* add 'elt' to 'lst' if it is not already in *)
let add_if_not_mem elt lst =
  if L.mem ~equal:(=) lst elt then
    lst
  else
    elt :: lst

(* remove all occurences of 'elt' from 'lst' *)
let remove_all elt lst =
  L.filter ~f:((<>) elt) lst

(* return the list of integers [start ; ... ; stop] *)
let list_of_ints start stop =
  let rec list_of_ints_priv start stop acc =
    if start = stop then
      stop :: acc
    else
      list_of_ints_priv start (stop - 1) (stop :: acc)
  in
  list_of_ints_priv start stop []

(* same as list_of_ints but convert them to floats afterwards *)
let list_of_floats start stop =
  L.map ~f:Float.of_int (list_of_ints start stop)

(* list of floats [start, start + step, ..., stop], without accumulation
   of floating point errors.
   step = (stop - start) / nb_steps *)
let discretize_range start nb_steps stop =
  let delta = stop -. start in
  let n     = itof nb_steps in
  let res   = ref []        in
  for i = nb_steps downto 0 do
    res := (((delta *. (itof i)) /. n) +. start) :: !res;
  done;
  !res

let shorter_list_first l1 l2 =
  compare (L.length l1) (L.length l2)

let longer_list_first l1 l2 =
  -(shorter_list_first l1 l2)

let string_of_list to_string sep l =
  "[" ^ S.concat ~sep:sep (L.map ~f:to_string l) ^ "]"

let rev_enrich_with f l =
  L.rev_map ~f:(fun x -> (f x, x)) l

exception Different_list_lengths;;

(* combine is not tail rec in the std lib *)
let combine lx ly =
  let rec combine_priv l1 l2 acc =
    match l1, l2 with
        x :: xs, y :: ys -> combine_priv xs ys ((x,y) :: acc)
      | []     , []      -> L.rev acc
      | _                -> raise Different_list_lengths
  in
  combine_priv lx ly []

let combine3 lx ly lz =
  let rec combine3_priv l1 l2 l3 acc =
    match l1, l2, l3 with
        x :: xs, y :: ys, z :: zs -> combine3_priv xs ys zs ((x,y,z) :: acc)
      | []     , []     , []      -> L.rev acc
      | _                         -> raise Different_list_lengths
  in
  combine3_priv lx ly lz []

exception List_length_not_multiple_of_three;;

let triplets_list_of_list l =
  let rec triplets_list_of_list_priv l acc =
    match l with
        x :: y :: z :: tail ->
          triplets_list_of_list_priv tail ((x,y,z) :: acc)
      | [] -> L.rev acc
      | _  -> raise List_length_not_multiple_of_three
  in
  triplets_list_of_list_priv l []

(* list iteration with a list of functions *)
let meta_iter ~funs l =
  L.iter
    ~f:(fun f -> L.iter ~f:f l)
    funs

(* product of all elements of a float list *)
let fprod l =
  L.fold l ~init:1.0 ~f:( *. )

(* product of all elements of an int list *)
let prod l =
  L.fold l ~init:1 ~f:( * )

let average_l l =
  let rec loop lst (sum, count) =
    match lst with
    | []      -> sum /. count
    | x :: xs -> loop xs (sum +. x, count +. 1.)
  in
  loop l (0., 0.)

(* standard deviation (sigma) *)
let std_dev l =
  let avg = average_l l in
  sqrt (average_l (L.map ~f:(fun x -> (x -. avg) *. (x -. avg)) l))

(* apply f to elements of l that satisfy p *)
let iter_on_some_elements p f l =
  L.iter
    ~f:(fun x -> if p x then f x else ())
    l

(* return a list of list by using a partitioning function *)
let group_by partition_f lst =
  let rec loop l acc = match l with
      [] -> L.rev acc
    | _  ->
        let curr, others = partition_f l in
        loop others (curr :: acc)
  in
  loop lst []

(* return Some (n first elements of lst) or None *)
let take_n n lst =
  let rec loop n l acc = match n with
      0 -> Some (L.rev acc)
    | _ -> match l with
          x :: xs -> loop (n - 1) xs (x :: acc)
        | _       -> None
  in
  loop n lst []

let take_with_rest n l =
  (BL.take n l, BL.drop n l)

let split_by n lst =
  if n < 1 then
    failwith (sprintf "my_utils.ml: split_by: invalid argument: %d" n)
  else
    let rec loop acc = function
      | [] -> BL.rev acc
      | l -> let n_elts, rest = take_with_rest n l in
             loop (n_elts :: acc) rest
    in loop [] lst

(* ============================== Sets ================================ *)

module IntSet = Int.Set

(* ============================== Time ================================ *)

let gettimeofday_fractional_part () =
  let tod     = Unix.gettimeofday() in
  let int_tod = Float.to_int tod    in
  tod -. (itof int_tod)

(* ============================== Tuples ============================== *)

(* compare 2 pairs based on their 1st element *)
let compare_fst (i1, _) (i2, _) =
  compare i1 i2

(* compare 2 pairs based on their 2nd element *)
let compare_snd (_, i1) (_, i2) =
  compare i1 i2

(* swap a pair *)
let swap (i, j) = (j, i)

(* create a sorted triple from a, b and c *)
let sort_triple a b c =
  match L.sort compare [a; b; c] with
    one :: two :: three :: [] -> (one, two, three)
  | _                         ->
      failwith "my_utils.ml: sort_triple: the impossible happened"

(* ============================== Hash tables ============================== *)

(* return all the values of a given hash table, unsorted *)
let values hash_table =
  HT.fold
    (fun _k v acc -> v :: acc)
    hash_table
    []

(* return all the keys of a given hash table, unsorted *)
let keys hash_table =
  HT.fold
    (fun k _v acc -> k :: acc)
    hash_table
    []

(* return all the (key, value) pairs of a given hash table, unsorted *)
let key_values hash_table =
  HT.fold
    (fun k v acc -> (k,v) :: acc)
    hash_table
    []

(* add key->value mapping to 'hash_table' if not already present *)
let add_if_not_present hash_table (key, value) =
  if not (HT.mem hash_table key) then
    HT.add hash_table key value

let hashtbl_of key_value_pairs =
  let res = HT.create (L.length key_value_pairs) in
  let add (k, v) = HT.add res k v in
  L.iter ~f:add key_value_pairs;
  res

(* transform a Hashtbl of key->value into a Hashtbl of value->key *)
let reverse_binding ht =
  let res = HT.create (HT.length ht) in
  HT.iter (fun key value -> HT.add res value key) ht;
  res

(* hash table lookup with key being printed out in exception if not found *)
let find key_printer ht key =
  let res = try HT.find ht key with
      Not_found -> failwith ("my_utils.ml: find: HT lookup failed for: "
                             ^ (key_printer key)) in
  res

(* set intersection on hash table keys *)
let inter ht1 ht2 =
  let l1 = HT.length ht1 in
  let l2 = HT.length ht2 in
  (* we'll iter on the smallest one *)
  let len, h1, h2 =
    if l1 < l2
    then l1, ht1, ht2
    else l2, ht2, ht1
  in
  let res = HT.create len in
  HT.iter
    (fun k v -> if HT.mem h2 k then HT.add res k v)
    h1;
  res

(* ============================== Maths ============================= *)

(* No Math.pi in Ocaml ?! *)
let pi       = 4.0 *. atan 1.0
let pi_div_2 = 2.0 *. atan 1.0

let is_NaN = Float.is_nan

(* flipping parameters *)
let flip f x y = f y x

let fabs x = if x < 0.0 then -.x
                        else   x

(* x =? v +/- epsilon *)
let around epsilon v x =
  (x >= v -. epsilon) &&
  (x <= v +. epsilon)

(* a famous function that does nothing *)
let identity x = x

let square x = x *. x

let cube x = x *. x *. x

let median_a xs =
  let n = A.length xs in
  A.sort xs ~compare:Float.compare;
  if n mod 2 = 1
  then xs.(n/2)
  else (xs.(n/2) +. xs.(n/2 - 1)) /. 2.0

let median_l xs =
  median_a (A.of_list xs)

(* significance of the correlation score c for a population of size n
   formula (14.6.2) for t at page 749 of numerical recipes 3rd ed. *)
let significance c n =
  let nf = Float.of_int n in
  c *. sqrt (( nf -. 2.) /. (1. -. c*.c))

(* approximation of the complementary error function erfc(x),
   comes from the book "numerical recipes" 3rd edition *)
let erfcc x =
  let   z = fabs x          in
  let   t = 2. /. (2. +. z) in
  let ans = t *. exp
    (-.z *. z -. 1.26551223 +. t *.
              (  1.00002368 +. t *.
              (  0.37409196 +. t *.
              (  0.09678418 +. t *.
              (-.0.18628806 +. t *.
              (  0.27886807 +. t *.
              (-.1.13520398 +. t *.
              (  1.48851587 +. t *.
              (-.0.82215223 +. t *.
                 0.17087277))))))))) in
  if x >= 0.0 then ans
              else 2.0 -. ans

(* Pearson correlation coefficient for float arrays, cross validated with some
   Python implementation of it that I have *)
let pearson_a a1 a2 =
  let tiny = 1.0e-20     in
  let    n = A.length a1 in
  let   nf = itof n      in
  if n <> A.length a2 then
    failwith "my_utils.ml: pearson_a: arrays length differ"
  else
    let p      = average_a a1 in
    let q      = average_a a2 in
    let sum_xx = ref 0.       in
    let sum_yy = ref 0.       in
    let sum_xy = ref 0.       in
    let process x' y' =
      let x    = x' -. p in
      let y    = y' -. q in
      let xx   = x *. x  in
      let yy   = y *. y  in
      let xy   = x *. y  in
        sum_xx := !sum_xx +. xx;
        sum_yy := !sum_yy +. yy;
        sum_xy := !sum_xy +. xy;
    in
    for i = 0 to n - 1 do
      process a1.(i) a2.(i);
    done;
    let r = !sum_xy /. (sqrt(!sum_xx *. !sum_yy) +. tiny)        in
    let z = 0.5 *. log((1.0 +. r +. tiny) /. (1.0 -. r +. tiny)) in
    (* approximation of Student's t probability valid for large n *)
    let t = erfcc(fabs(z *. sqrt(nf -. 1.0)) /. 1.4142136)       in
    (r, t)

(* comes from Biocaml, not me *)
let spearman_rank arr =
  let arr = A.copy arr in
  let arr = A.mapi arr (fun i a -> a,i) in
  A.sort arr ~compare:(fun (a,_) (b,_) -> Pervasives.compare a b);
  let g _prev il ans =
    let count = L.length il in
    let n = count + (L.length ans) in
    let hi = Float.of_int n in
    let lo = Float.of_int (n - count + 1) in
    let rank = (hi +. lo) /. 2. in
    L.append (L.map ~f:(fun i -> rank,i) il) ans
  in
  let f (prev, il, ans) (x,i) =
    let count = L.length il in
    if count = 0 then
      x, [i], ans
    else if x = prev then
      x, i::il, ans
    else
      x, [i], g prev il ans
  in
  let prev, il, ans = fold_left ~f:f ~init:(0.,[],[]) arr in
  let ans = g prev il ans in
  let ans = L.sort (fun (_,a) (_,b) -> Pervasives.compare a b) ans in
  A.of_list (L.map ~f:fst ans)

(* Spearman rank-order correlation coefficient *)
let spearman_a (a1:float array) (a2:float array) =
  pearson_a (spearman_rank a1) (spearman_rank a2)

let spearman_l l1 l2 =
  spearman_a (A.of_list l1) (A.of_list l2)

let pearson_l l1 l2 =
  pearson_a (A.of_list l1) (A.of_list l2)

(* initialize the RNG using system time in ms as the seed (returned) *)
let init_RNG () =
  let ms_int = ftoi (1000. *. gettimeofday_fractional_part()) in
  RNG.init ms_int;
  ms_int

(* ============================= Options ============================ *)

(* the unsafe get *)
let get_opt_or_fail ~opt ~msg = match opt with
    Some v -> v
  | None   -> failwith msg

(* ============================= System / UNIX ============================ *)

let detect_cmd cmd =
  match Sys.command ("which " ^ cmd ^ " > /dev/null") with
      0 -> true
    | _ -> false

(* fail if command is not available *)
let require_cmd cmd =
  match Sys.command ("which " ^ cmd ^ " > /dev/null") with
      0 -> () (* command is here *)
    | _ -> failwith ("my_utils.ml: require_cmd: command not found: " ^ cmd)

let getenv_or_fail variable_name failure_message =
  get_opt_or_fail ~opt:(Sys.getenv variable_name) ~msg:failure_message

exception Command_failed of string;;

(* run a command,
   ~debug  -> dry-run only
   ~strict -> raise (Command_failed msg) on non 0 exit status of the command *)
let run_command ?(debug = false) ?(strict = true) cmd =
  P.printf "running:\n%s\n%!" cmd;
  if not debug then
    let ret_code = Sys.command cmd in
    (* P.printf "ret_code: %d\n%!" ret_code; *)
    match ret_code with
        0 -> ()
      | _ -> if strict then raise (Command_failed ("command failed: " ^ cmd))

(* gzip a file using the given compression level
   return the compressed filename
   the original file is removed *)
let gzip filename comp_level =
  run_command (P.sprintf "gzip -f -%d %s" comp_level filename);
  (filename ^ ".gz")

(* gunzip a file
   return the uncompressed filename
   the original file can be kept but is removed by default *)
let gunzip ?(keep = false) filename =
  let res = F.chop_extension filename in
  let cmd =
    if keep
    then P.sprintf "gunzip -c %s > %s" filename res
    else P.sprintf "gunzip -f %s"      filename
  in
  run_command cmd;
  res

(* run the given command and return its output as a string *)
let get_command_output cmd =
  P.printf "running:\n%s\n%!" cmd;
  let cmd_out   = Unix.open_process_in cmd in
  let buff_size = 1024*1024                in
  let buff      = Buffer.create buff_size  in
  let _         =
    try
      while true do
        Buffer.add_string buff (input_line cmd_out);
        Buffer.add_char   buff '\n'; (* was stripped by input_line *)
      done;
    with
    | End_of_file -> ignore(Unix.close_process_in cmd_out)
    | exn         -> ignore(Unix.close_process_in cmd_out); raise exn
  in
  Buffer.contents buff

(* run the given command and return its output as a list of lines *)
let get_command_output_lines cmd =
  P.printf "running:\n%s\n%!" cmd;
  let cmd_out = Unix.open_process_in cmd in
  let buff    = ref [] in
  let _       =
    try
      while true do
        buff := (input_line cmd_out) :: !buff;
      done;
    with
    | End_of_file -> ignore(Unix.close_process_in cmd_out)
    | exn         -> ignore(Unix.close_process_in cmd_out); raise exn
  in
  L.rev !buff

(* ============================== Files ============================== *)

(* get_extension "toto.txt" -> "txt"
   get_extension "toto"     -> ""
   get_extension ""         -> ""    *)
let get_extension filename =
  match F.split_extension filename with
      _, None     -> ""
    | _, Some ext -> ext

(* check filename has one in a list of possible extensions *)
let enforce_any_file_extension filename possible_exts_list =
  let res =
    L.filter
      ~f:(fun ext -> F.check_suffix filename ext)
      possible_exts_list
  in
  match res with
    [] -> failwith ("my_utils.ml: enforce_any_file_extension: " ^
                    filename ^ " is not any of " ^
                    (string_of_list identity "; " possible_exts_list))
  | _  -> ()

let filename_with_different_extension filename previous_ext new_ext =
  (* don't chop silently if the extension is not the expected one,
     crash loud and early instead *)
  enforce_any_file_extension filename [previous_ext];
  (F.chop_suffix filename previous_ext) ^ new_ext

(* test if a valid gzip file extension is present in filename *)
let is_a_compressed_file filename =
  F.check_suffix filename ".gz"   ||
  F.check_suffix filename ".gzip" ||
  F.check_suffix filename ".Z"

(* marshal a value to a possibly compressed file, determined based on the
   extension of f *)
let marshal v f =
  let marshal_priv v f =
    let out = open_out f in
    Mar.to_channel out v [];
    close_out_noerr out
  in
  if is_a_compressed_file f then begin
    let uncompressed_f = F.chop_extension f in
    marshal_priv v uncompressed_f;
    ignore(gzip uncompressed_f 1);
  end else
    marshal_priv v f

(* unmarshal a value from a possibly compressed file *)
let unmarshal f =
  let unmarshal_priv f =
    let input = open_in f              in
    let v     = Mar.from_channel input in
    close_in_noerr input;
    v
  in
  if is_a_compressed_file f then
    let uncompressed_f = gunzip ~keep:true f           in
    let res            = unmarshal_priv uncompressed_f in
    Sys.remove uncompressed_f;
    res
  else
    unmarshal_priv f

(* ============================= My error monad ======================== *)

(* the very badly named "return" in Haskell, it is just a constructor in fact *)
type ('a, 'b) result =
  | Result of 'a
  | Error  of 'b

(* "bind" in Haskell *)
let try_to_apply m f = match m with
  | Result r -> f r (* apply the function *)
  | Error  _ -> m   (* leave the value untouched and propagate it downstream *)

(* "bind" as an operator *)
let (>>=) m f = try_to_apply m f
