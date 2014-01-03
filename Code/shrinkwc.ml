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

let nonzero_sign a = if a > 0 then 1 else -1;;

let output_matrix matrix file_name =
  for k = 1 to Array.length matrix do
    for j = 1 to Array.length matrix.(0) do
      output_byte file_name (int_of_float matrix.(k-1).(j-1))
    done;
  done;
;;

let rec wc2matrix matrix buffer info counter =
  if counter = info.num_of_codes then
    matrix
  else
    let code = Scanf.bscanf buffer "%i %i %i %i %i\n" retrieve_code in
    let coefficient = info.threshold *. (2.0 ** ((float (2 * code.level)) /.  info.p)) *. (float code.value) -. (float (nonzero_sign code.value)) *. info.threshold in
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
  wc2matrix matrix buffer info (counter+1);;

let coeffs_id = Sys.argv.(1);;
let buffer    = Scanf.Scanning.from_file coeffs_id;;
let info_line = Scanf.bscanf buffer "%s %i %f %f %i\n" retrieve_info;;

let dimension = int_of_float (2.0 ** (float info_line.dyadic_level));;
let dim_str   = string_of_int dimension;;
let matrix    = Array.make_matrix dimension dimension 0.0;;
let file_name = coeffs_id^".pgm";;
let file      = open_out_bin file_name;;
let header    = "P5\n"^dim_str^" "^dim_str^"\n255\n";;

output_string file header;;
output_matrix (wc2matrix matrix buffer info_line 0) file;;
close_out file;;
