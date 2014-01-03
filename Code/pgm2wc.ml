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

(* MAIN FUNCTION *)

let image_id     = Sys.argv.(1);;
let p            = float_of_string Sys.argv.(2);;
let threshold    = float_of_string Sys.argv.(3);;
let coeffs_file  = Sys.argv.(4);;
let image        = open_in_bin image_id;;
let magic_number = input_line image;;
let dimensions   = input_line image;;
let grey_levels  = input_line image;;
(* If grey_levels >= 256, then we read two bytes instead of one.  
   I'll do that later. *)

  let step1        = String.index dimensions ' ';;
  let image_length = int_of_string (String.sub dimensions 0 step1);;
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

Unix.system "/usr/bin/mkdir tmp";;
let temp_coeffs   = "tmp/"^coeffs_file;;
let coefficients  = open_out temp_coeffs;;
compute_wavelet_coefficients img 
			     image_length 
			     level 
			     p
			     threshold
			     coefficients;;
(* And before we close "coefficients", we must include in the info line the
number of lines (= number of coefficients included) *)

close_out coefficients;;

let info_line     = open_out "tmp/info";;
Printf.fprintf info_line "%s %i %f %f " image_id level p threshold;;
close_out info_line;;

let command1 = "wc -l tmp/"^coeffs_file^" | awk '{ print $1 }' >> tmp/info";;
let command2 = "cat tmp/info tmp/"^coeffs_file^" > "^coeffs_file;;
Unix.system command1;;
Unix.system command2;;
Unix.system "/usr/bin/rm -r tmp";;
