(*
				     pgm2wc 

									      *)

let send_code coeff int1 int2 dyadic_level wavelet_indicator p threshold file =
  let weight = (2.0 ** (2.0 *. (float dyadic_level) /. p)) *. threshold in
  let code   = int_of_float(coeff /. weight) in
  if code != 0 then
    Printf.fprintf file "%i %i %i %i %i\n" code int1 int2 dyadic_level wavelet_indicator;;

let rec compute_wavelet_coefficients matrix 
				     matrix_length 
				     dyadic_level 
				     p
				     threshold
				     file =
  if dyadic_level = 0 then 
    send_code matrix.(0).(0) 0 0 0 0 p threshold file
  else begin
    let proj_coeffs = Array.make_matrix (matrix_length / 2) (matrix_length / 2) 0.0 in
    for k = 1 to (matrix_length / 2) do
      for j = 1 to (matrix_length / 2) do
		let weigth = (2.0 ** (2.0 *. (float dyadic_level) /. p)) *. threshold in
		let a = matrix.(2*(k-1)).(2*(j-1)) in
		let b = matrix.(2*(k-1)).(2*(j-1)+1) in
		let c = matrix.(2*(k-1)+1).(2*(j-1)) in
		let d = matrix.(2*(k-1)+1).(2*(j-1)+1) in
		let c1 = (a+.c-.b-.d) /. 4.0 in
		let c2 = (c+.d-.a-.b) /. 4.0 in
		let c3 = (b+.c-.a-.d) /. 4.0 in
		let c4 = (a+.b+.c+.d) /. 4.0 in
		send_code c1 (k-1) (j-1) (dyadic_level-1) 1 p threshold file;
		send_code c2 (k-1) (j-1) (dyadic_level-1) 2 p threshold file;
		send_code c3 (k-1) (j-1) (dyadic_level-1) 3 p threshold file;
		proj_coeffs.(k-1).(j-1) <- c4
	done;
      done;
      compute_wavelet_coefficients proj_coeffs 
				   (matrix_length / 2) 
				   (dyadic_level - 1) 
				   p
				   threshold
				   file
  end;;
let pgm2wc image_id img image_length level p threshold wc_file =
  Unix.system "/usr/bin/mkdir tmp";
  let temp_coeffs = "tmp/"^wc_file in
  let coefficients = open_out temp_coeffs in
  compute_wavelet_coefficients img 
			       image_length 
			       level 
			       p
			       threshold
			       coefficients;
  close_out coefficients;
  let info_file = open_out "tmp/info" in
  Printf.fprintf info_file "%s %i %f %f " image_id level p threshold;
  close_out info_file;
  let command1 = "wc -l tmp/"^wc_file^" | awk '{ print $1 }' >> tmp/info" in
  let command2 = "cat tmp/info tmp/"^wc_file^" > "^wc_file in
  Unix.system command1;
  Unix.system command2;
  Unix.system "/usr/bin/rm -r tmp";;


(* Utilities to deal with wc files *)
let retrieve_number_of_coefficients a b c d e = e;;

let read_number_of_coefficients wc_file =
  let buffer = Scanf.Scanning.from_file wc_file in
  Scanf.bscanf buffer "%s %f %f %f %f" retrieve_number_of_coefficients;;


(*

			       wc2pgm & wc2shrink
									      *)

type info = 
{
  image_id:     string;
  dyadic_level: int;
  p:            float;
  threshold:    float;
  num_of_codes: int
};;

type code =
{
  value:   int;
  index1:  int;
  index2:  int;
  level:   int;
  wavelet: int;
}

let retrieve_info s i1 f1 f2 i2 =
{
  image_id =     s;
  dyadic_level = i1;
  p =            f1;
  threshold =    f2;
  num_of_codes = i2
};;

let retrieve_code i1 i2 i3 i4 i5 =
{
  value =   i1;
  index1 =  i2;
  index2 =  i3;
  level =   i4;
  wavelet = i5
};;

let int2float_sign a = if a > 0 then 1.0 else -1.0;;

let output_matrix matrix file_name =
  for k = 1 to Array.length matrix do
    for j = 1 to Array.length matrix.(0) do
      output_byte file_name (int_of_float matrix.(k-1).(j-1))
    done;
  done;
;;

let rec wc2matrix matrix buffer info shrinkage_coeff counter =
  if counter = info.num_of_codes then
    matrix
  else
    let code = Scanf.bscanf buffer "%i %i %i %i %i\n" retrieve_code in
    let coefficient = info.threshold *. (2.0 ** ((float (2 * code.level)) /.  info.p)) *. (float code.value) -. (float shrinkage_coeff) *.  (int2float_sign code.value) *. info.threshold in
    let depth = int_of_float (2.0 ** (float (info.dyadic_level - code.level))) in
    if code.wavelet == 0 then
      begin
	for k = 1 to depth do
	  for j = 1 to depth do
	    let value = matrix.(k-1).(j-1) in
	    matrix.(k-1).(j-1) <- value +. coefficient
	  done;
	done
      end;
    if code.wavelet == 1 then
      begin
	for k = 1 to depth do
	  for j = 1 to (depth/2) do
	    let value = matrix.(depth*code.index1+k-1).(depth*code.index2+j-1) in
	    matrix.(depth*code.index1+k-1).(depth*code.index2+j-1) <- value +.
	    coefficient
	  done;
	done;
	for k = 1 to depth do
	  for j = (depth/2)+1 to depth do
	    let value = matrix.(depth*code.index1+k-1).(depth*code.index2+j-1) in
	    matrix.(depth*code.index1+k-1).(depth*code.index2+j-1) <- value -. coefficient
	  done;
	done;
      end;
    if code.wavelet == 2 then
      begin
	for k = 1 to (depth/2) do
	  for j = 1 to depth do
	    let value = matrix.(depth*code.index1+k-1).(depth*code.index2+j-1) in
	    matrix.(depth*code.index1+k-1).(depth*code.index2+j-1) <- value -.
	    coefficient
	  done;
	done;
	for k = (depth/2)+1 to depth do
	  for j = 1 to depth do
	    let value = matrix.(depth*code.index1+k-1).(depth*code.index2+j-1) in
	    matrix.(depth*code.index1+k-1).(depth*code.index2+j-1) <- value +. coefficient
	  done;
	done;
      end;
    if code.wavelet == 3 then
      begin
	for k = 1 to (depth/2) do
	  for j = 1 to (depth/2) do
	    let value = matrix.(depth*code.index1+k-1).(depth*code.index2+j-1) in
	    matrix.(depth*code.index1+k-1).(depth*code.index2+j-1) <- value -.
	    coefficient
	  done;
	done;
	for k = (depth/2)+1 to depth do
	  for j = 1 to (depth/2) do
	    let value = matrix.(depth*code.index1+k-1).(depth*code.index2+j-1) in
	    matrix.(depth*code.index1+k-1).(depth*code.index2+j-1) <- value +. coefficient
	  done;
	done;
	for k = 1 to (depth/2) do
	  for j = (depth/2)+1 to depth do
	    let value = matrix.(depth*code.index1+k-1).(depth*code.index2+j-1) in
	    matrix.(depth*code.index1+k-1).(depth*code.index2+j-1) <- value +.
	    coefficient
	  done;
	done;
	for k = (depth/2)+1 to depth do
	  for j = (depth/2)+1 to depth do
	    let value = matrix.(depth*code.index1+k-1).(depth*code.index2+j-1) in
	    matrix.(depth*code.index1+k-1).(depth*code.index2+j-1) <- value -. coefficient
	  done;
	done;
      end; 
  wc2matrix matrix buffer info shrinkage_coeff (counter+1);;

let wc2pgm coeffs_id =
  let buffer    = Scanf.Scanning.from_file coeffs_id in
  let info_line = Scanf.bscanf buffer "%s %i %f %f %i\n" retrieve_info in

  let dimension = int_of_float (2.0 ** (float info_line.dyadic_level)) in
  let dim_str   = string_of_int dimension in
  let matrix    = Array.make_matrix dimension dimension 0.0 in
  let file_name = coeffs_id^".pgm" in
  let file      = open_out_bin file_name in
  let header    = "P5\n"^dim_str^" "^dim_str^"\n255\n" in

  output_string file header;
  output_matrix (wc2matrix matrix buffer info_line 0 0) file;
  close_out file;;

let wc2shrink coeffs_id =
  let buffer    = Scanf.Scanning.from_file coeffs_id in
  let info_line = Scanf.bscanf buffer "%s %i %f %f %i\n" retrieve_info in

  let dimension = int_of_float (2.0 ** (float info_line.dyadic_level)) in
  let dim_str   = string_of_int dimension in
  let matrix    = Array.make_matrix dimension dimension 0.0 in
  let file_name = coeffs_id^".pgm" in
  let file      = open_out_bin file_name in
  let header    = "P5\n"^dim_str^" "^dim_str^"\n255\n" in

  output_string file header;
  output_matrix (wc2matrix matrix buffer info_line 1 0) file;
  close_out file

(*

	       L_2 norm of the difference between two pgm images
									      *)
let nimi image_id1 image_id2 =
  let image1       = open_in_bin image_id1 in
  let image2       = open_in_bin image_id2 in
  let magic_number = input_line image1 in
  input_line image2;
  let dimensions   = input_line image1 in
  input_line image2;
  let grey_levels  = input_line image1 in
  input_line image2;

    let step         = String.index dimensions ' ' in
    let image_length = int_of_string (String.sub dimensions 0 step) in

  let norm = [|0.0|] in
  for k = 1 to image_length do
    for j = 1 to image_length do
      let value = ((input_byte image1) - (input_byte image2)) in
      let previous = norm.(0) in
      norm.(0) <- previous +. (float (value * value)) /. (float (image_length * image_length))
    done;
  done;
  sqrt norm.(0);;

(*
				 Balphatau_norm
									      *)
let pi     = 3.14159265358979;;

let float_abs a = if a < 0.0 then (a *. -1.0) else a;;

let balphatau_norm_helper file tau =
  let norm = Array.make 1 0.0 in
  let buffer = Scanf.Scanning.from_file file in
  let info_line = Scanf.bscanf buffer "%s %i %f %f %i\n" retrieve_info in
  for k = 1 to info_line.num_of_codes do
    let code = Scanf.bscanf buffer "%i %i %i %i %i\n" retrieve_code in
    let coefficient = info_line.threshold *. (2.0 ** ((float (2 * code.level)) /. info_line.p)) *. (float code.value) in
    let value = norm.(0) in
    norm.(0) <- value +. ((float_abs coefficient) ** tau)
  done;
  norm.(0) ** (1.0 /. tau);;

(*
			       Least Squares Line
									      *)
let identity a = a;;
let square a = a *. a;;

let rec sum_filtered_column matrix column filter counter size value =
  if counter > size then 
    value 
  else
    let new_value = (filter matrix.(counter-1).(column)) in
    sum_filtered_column matrix column filter (counter+1) size (value+.new_value);;

let rec sum_product_of_columns matrix counter size value =
  if counter > size then
    value
  else 
    let new_value = matrix.(counter-1).(0) *. matrix.(counter-1).(1) in
    sum_product_of_columns matrix (counter+1) size (value+.new_value);;


(*   MAIN   *)
let image_id     = Sys.argv.(1);;
let image        = open_in_bin image_id;;
let magic_number = input_line image;;
let dimensions   = input_line image;;
let grey_levels  = input_line image;;
(* If grey_levels >= 256, then we read two bytes instead of one.  
   I'll do that later. *)

  let step         = String.index dimensions ' ';;
  let image_length = int_of_string (String.sub dimensions 0 step);;
  let level        = int_of_float (log (float image_length) /. (log 2.0));;
  (* We assume that the input images are square.  
     I'll do the general case later too. *)

(* Read the image into a matrix *)
let img = Array.make_matrix image_length image_length 0.0;;
for k = 1 to image_length do
  for j = 1 to image_length do
    img.(k-1).(j-1) <- float_of_int (input_byte image)
  done;
done;;

(*   This computes an approximation to \alpha  *)
let data = open_out "data";;

(* Create the pairs (-log N, log ||f-fN||_2)   *)
let pairs = Array.make_matrix 25 2 0.0;;

for k = 1 to 25 do
  let wc_file   = "wco"^(string_of_int k) in
  let pgm_file  = wc_file^".pgm" in
  let threshold = (float (k-1)) /. 12000.0 +. 0.0001 in
  pgm2wc image_id img image_length level 2.0 threshold wc_file;
  wc2pgm wc_file;
  let minus_log_num_coefs = -.log (read_number_of_coefficients wc_file) in
  let log_error           =   log (nimi image_id pgm_file) in
  pairs.(k-1).(0) <- minus_log_num_coefs;
  pairs.(k-1).(1) <- log_error;
  Printf.fprintf data "%f %f\n" minus_log_num_coefs log_error
done;;
close_out data;;

let sum_x  = (sum_filtered_column pairs 0 identity 1 25 0.0);;
let sum_y  = (sum_filtered_column pairs 1 identity 1 25 0.0);;
let sum_x2 = (sum_filtered_column pairs 0 square   1 25 0.0);;
let sum_y2 = (sum_filtered_column pairs 1 square   1 25 0.0);;
let sum_xy = (sum_product_of_columns pairs 1 25 0.0);;

let denom = 25.0 *. sum_x2 -. (square sum_x);;
let slope = (25.0 *. sum_xy -. sum_x *. sum_y) /. denom;; 
let b     = (sum_y *. sum_x2 -. sum_x *. sum_xy) /. denom;;
let alpha = slope /. 2.0;;
let tau   = 2.0 /. (1.0 +. alpha);;
let balphatau_norm = balphatau_norm_helper "wco1" tau;;
(* *)
  print_string image_id;;
  print_string "\nalpha: ";;
  print_float alpha;;
  print_string "\ntau:   ";;
  print_float tau;;
  print_string "\nBat:   ";;
  print_float balphatau_norm;;
(* *)

Unix.system "echo 'set terminal postscript' > script";;
let gnuplot_command = "echo 'plot \"data\", "^(string_of_float b)^"+"^(string_of_float slope)^"*x' >> script";;
Unix.system gnuplot_command;;
Unix.system "gnuplot script > plot.eps";;

(* We add noise N(0,32) into the image-matrix, 
   and perform some de-noising on it:  *)

let noise_buffer = open_out_bin "noise.pgm";;
let stdev = 16.0;;
Random.self_init;;

for k = 1 to image_length do
  for j = 1 to image_length do
    let random_number = Random.float 1.0 in (* pseudo-random generated. *)
    (* the noise is normal distributed N(0,1), and then multiplied by stdev. *)
    let noise = stdev *. (sqrt (-2.0 *. (log random_number))) *. (cos (2.0 *. pi *. random_number)) in
    let value = img.(k-1).(j-1) in
    img.(k-1).(j-1) <- max 0.0 (min 255.0 (value +. noise))
  done;
done;;

output_string noise_buffer "P5\n";
output_string noise_buffer (string_of_int image_length);
output_string noise_buffer " ";
output_string noise_buffer (string_of_int image_length);
output_string noise_buffer "\n255\n";
for k = 1 to image_length do
  for j = 1 to image_length do
    output_byte noise_buffer (int_of_float img.(k-1).(j-1))
  done;
done;;
close_out noise_buffer;;

(* de-noising step: we perform VisuShrink and EasyShrink *)

let m = (float (image_length * image_length));;
let visushrink = stdev *. sqrt ( 2.0 *. (log m))  /. (float image_length);;
let easyshrink = stdev *. (sqrt (float_abs ((2.0 -. tau) *. (log m) -. 2.0 *. tau *. (log (balphatau_norm /. stdev))))) /. (float image_length);;

(* *)
  print_string "\nVisu:  ";;
  print_float visushrink;;
  print_string "\nEasy:  ";;
  print_float easyshrink;;
(* *)

let visushrink_wc_file = "VisuShrink";;
let easyshrink_wc_file = "EasyShrink";;
pgm2wc "noise.pgm" img image_length level 2.0 visushrink visushrink_wc_file;
wc2shrink visushrink_wc_file;;
pgm2wc "noise.pgm" img image_length level 2.0 easyshrink easyshrink_wc_file;
wc2shrink easyshrink_wc_file;;

(* *)
  print_string "\nerror VisuShrink: ";;
  print_float (nimi image_id "VisuShrink.pgm");;
  print_string "\nerror EasyShrink: ";;
  print_float (nimi image_id "EasyShrink.pgm");;
  print_string "\n";;
(* *)
