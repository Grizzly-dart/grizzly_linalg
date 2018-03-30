import 'package:grizzly_array/grizzly_array.dart';
import 'package:grizzly_linalg/grizzly_linalg.dart';
import 'package:test/test.dart';

void main() {
  group('A group of tests', () {
    test('Factorize', () {
      final a = new Double2D([
        [1.0, 4.0],
        [2.0, 5.0],
      ]);
      final LU xlu = lu(a);

      expect(
          xlu.p,
          new Int2D([
            [0, 1],
            [1, 0]
          ]));

      expect(
          xlu.l,
          new Double2D([
            [1.0, 0.0],
            [0.5, 1.0],
          ]));

      expect(
          xlu.u,
          new Double2D([
            [2.0, 5.0],
            [0.0, 1.5]
          ]));

      expect(xlu.p.matmul(xlu.l).matmul(xlu.u), a);
    });

    test('Solve', () {
      final a = new Double2D([
        [1.0, 4.0],
        [2.0, 5.0],
      ]);
      final b = a.matmul(new Double2D.aCol(<int>[5, 2]));
      final LU xlu = lu(a);

      expect(
          xlu.solve(b),
          ints2([
            [5],
            [2]
          ]));
    });
  });
}
