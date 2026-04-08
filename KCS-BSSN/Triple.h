#pragma once
template <class T1, class T2, class T3>
class Triple {
public:
	T1 first;
	T2 second;
	T3 third;
	Triple(T1 f, T2 s, T3 t) {
		first = f;
		second = s;
		third = t;
	}
};