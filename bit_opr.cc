//NOTE: numtype should be defined before including this file
// a is always the index of the greatest bit
// b is always the index of the smallest bit
template<int b>
struct bit{enum{set=numtype(1)<<b,unset = !set};};
template<int a, int b>
struct bits_in{enum{set=(numtype(1)<<a)|bits_in<a-1,b>::set,unset = !set};};
template<int b>
struct bits_in<b,b>{enum{set=bit<b>::set,unset = !set};};
template<int a, int b>
inline numtype occupy_bits(numtype small_value){
	numtype cut=small_value&bits_in<a-b,0>::set;
	if(small_value!=cut)throw;
	return cut<<b;
}
