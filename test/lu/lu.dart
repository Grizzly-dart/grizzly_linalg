import 'package:grizzly_series/grizzly_series.dart';
import 'package:grizzly_linalg/grizzly_linalg.dart';
import 'package:test/test.dart';

void main() {
	group('A group of tests', () {
		setUp(() {
		});

		test('First Test', () {
			final x = array2D([
				[1, 2],
				[2, 3],
				[3, 4],
				[4, 5],
				[5, 6],
			]);
			final y = (x * [5, 2]).row.sum + 7;
			print(y);

			final LU xLU = lu(x);

			print(xLU.pivotMatrix);

			print(xLU.lowerFactor);

			print(xLU.upperFactor);

			print(xLU.pivotMatrix * xLU.lowerFactor * xLU.upperFactor);
		});
	});
}