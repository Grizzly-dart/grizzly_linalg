library grizzly.linalg.norm;

import 'dart:math' as math;

double norm1(Iterable<num> v) {
	if(v.length == 0) return 0.0;
	if(v.length == 1) return v.first.abs().toDouble();

	double sum = 0.0;
	for(num s in v) sum += s.abs();
	return sum;
}

double norm2(Iterable<num> v) {
	if(v.length == 0) return 0.0;
	if(v.length == 1) return v.first.abs().toDouble();

	// TODO improve efficiency
	double sum = 0.0;
	for(num s in v) sum += s * s;
	return math.sqrt(sum);
}