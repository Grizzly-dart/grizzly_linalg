import 'package:grizzly_array/grizzly_array.dart';
import 'package:grizzly_linalg/grizzly_linalg.dart';
import 'package:test/test.dart';

void main() {
	group('A group of tests', () {
		setUp(() {
		});

		test('First Test', () {
			final x = double2D([
				[1.0, 2.0],
				[2.0, 3.0],
				[3.0, 4.0],
				[4.0, 5.0],
				[5.0, 6.0],
			]);
			final y = (x * [5, 2]).row.sum + 7;
			print(y.asIterable);

			final LU xlu = lu(x);

			print(xlu.pivotMatrix);

			print(xlu.lowerFactor);

			print(xlu.upperFactor);

			print(xlu.pivotMatrix * xlu.lowerFactor * xlu.upperFactor);
		});
	});
}