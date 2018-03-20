import 'package:grizzly_array/grizzly_array.dart';
import 'package:grizzly_linalg/grizzly_linalg.dart';
import 'package:test/test.dart';

void main() {
  group('A group of tests', () {
    setUp(() {});

    test('QR 0', () {
      final x = array2D([
        [1, 2],
        [2, 3],
        [3, 4],
        [4, 5],
        [5, 6],
      ]);
      final y = (x * [5, 2]).row.sum + 7;

      final xQR = qr(x);
      print(xQR.q);
      print(xQR.r);

      final b = solve(xQR.r, xQR.q.transpose * y.transpose);
      print(b);
    });

    test('QR 1', () {
      final x = array2D([
        [1],
        [2],
        [3],
        [4],
        [5],
      ]);

      final y = x.col[0] * 5;

      final xQR = qr(x);
      print(xQR.q);
      print(xQR.r);

      print(xQR.q.dot(xQR.r));

      print(xQR.solve(y.transpose));
    });
  });
}
