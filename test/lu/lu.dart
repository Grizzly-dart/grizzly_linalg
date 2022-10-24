import 'package:grizzly/grizzly.dart';
import 'package:grizzly_linalg/grizzly_linalg.dart';
import 'package:test/test.dart';

void main() {
  group('LU', () {
    group('Singularity', () {
      group('NonSingular', () {
        test('1', () {
          final a = [
            [1.0, 0.0],
            [0.0, 1.0]
          ];
          final LU xlu = lu(a);

          expect(xlu.p.matmul(xlu.l).matmul(xlu.u), a);
          expect(xlu.isSingular, false);
        });

        test('2', () {
          final a = [
            [0.0, 1.0],
            [1.0, 0.0],
          ];
          final LU xlu = lu(a);

          expect(xlu.p.matmul(xlu.l).matmul(xlu.u), a);
          expect(xlu.isSingular, false);
        });

        test('3', () {
          final a = [
            [1.0, 1.0],
            [0.0, 1.0]
          ];
          final LU xlu = lu(a);

          expect(xlu.p.matmul(xlu.l).matmul(xlu.u), a);
          expect(xlu.isSingular, false);
        });

        test('4', () {
          final a = [
            [0.0, 1.0],
            [1.0, 1.0]
          ];
          final LU xlu = lu(a);

          expect(xlu.p.matmul(xlu.l).matmul(xlu.u), a);
          expect(xlu.isSingular, false);
        });

        test('5', () {
          final a = [
            [1.0, 0.0],
            [1.0, 1.0]
          ];
          final LU xlu = lu(a);

          expect(xlu.p.matmul(xlu.l).matmul(xlu.u), a);
          expect(xlu.isSingular, false);
        });

        test('6', () {
          final a = [
            [1.0, 1.0],
            [1.0, 0.0]
          ];
          final LU xlu = lu(a);

          expect(xlu.p.matmul(xlu.l).matmul(xlu.u), a);
          expect(xlu.isSingular, false);
        });
      });

      group('Singular', () {
        test('0', () {
          final a = [
            [0.0, 0.0],
            [0.0, 0.0]
          ];
          final LU xlu = lu(a);

          expect(xlu.p.matmul(xlu.l).matmul(xlu.u), a);
          expect(xlu.isSingular, true);
        });

        test('1', () {
          final a = [
            [1.0, 0.0],
            [0.0, 0.0]
          ];
          final LU xlu = lu(a);

          expect(xlu.p.matmul(xlu.l).matmul(xlu.u), a);
          expect(xlu.isSingular, true);
        });

        test('2', () {
          final a = [
            [0.0, 1.0],
            [0.0, 0.0]
          ];
          final LU xlu = lu(a);

          expect(xlu.p.matmul(xlu.l).matmul(xlu.u), a);
          expect(xlu.isSingular, true);
        });

        test('3', () {
          final a = [
            [0.0, 0.0],
            [0.0, 1.0]
          ];
          final LU xlu = lu(a);

          expect(xlu.p.matmul(xlu.l).matmul(xlu.u), a);
          expect(xlu.isSingular, true);
        });

        test('4', () {
          final a = [
            [0.0, 0.0],
            [1.0, 0.0]
          ];
          final LU xlu = lu(a);

          expect(xlu.p.matmul(xlu.l).matmul(xlu.u), a);
          expect(xlu.isSingular, true);
        });

        test('5', () {
          final a = [
            [1.0, 1.0],
            [0.0, 0.0]
          ];
          final LU xlu = lu(a);

          expect(xlu.p.matmul(xlu.l).matmul(xlu.u), a);
          expect(xlu.isSingular, true);
        });

        test('6', () {
          final a = [
            [0.0, 1.0],
            [0.0, 1.0]
          ];
          final LU xlu = lu(a);

          expect(xlu.p.matmul(xlu.l).matmul(xlu.u), a);
          expect(xlu.isSingular, true);
        });

        test('7', () {
          final a = [
            [0.0, 0.0],
            [1.0, 1.0]
          ];
          final LU xlu = lu(a);

          expect(xlu.p.matmul(xlu.l).matmul(xlu.u), a);
          expect(xlu.isSingular, true);
        });

        test('8', () {
          final a = [
            [1.0, 0.0],
            [1.0, 0.0]
          ];
          final LU xlu = lu(a);

          expect(xlu.p.matmul(xlu.l).matmul(xlu.u), a);
          expect(xlu.isSingular, true);
        });

        test('9', () {
          final a = [
            [1.0, 1.0],
            [1.0, 1.0]
          ];
          final LU xlu = lu(a);

          expect(xlu.p.matmul(xlu.l).matmul(xlu.u), a);
          expect(xlu.isSingular, true);
        });
      });
    });

    group('Factorization', () {
      test('0', () {
        final a = [
          [1.0, 4.0],
          [2.0, 5.0],
        ];
        final LU xlu = lu(a);

        expect(xlu.p, [
          [0, 1],
          [1, 0]
        ]);

        expect(
            xlu.l,
            [
              [1.0, 0.0],
              [0.5, 1.0],
            ]);

        expect(
            xlu.u,
            [
              [2.0, 5.0],
              [0.0, 1.5]
            ]);

        expect(xlu.p.matmul(xlu.l).matmul(xlu.u), a);
        expect(xlu.isSingular, false);
      });

      test('1', () {
        final a = [
          [6.0, 1.0],
          [4.0, -2.0],
          [2.0, 8.0]
        ];
        final LU xlu = lu(a);

        expect(xlu.p.matmul(xlu.l).matmul(xlu.u), a);
        expect(() => xlu.isSingular, throwsUnsupportedError);
      });

      test('2', () {
        final a = [
          [4.0, 3.0],
          [6.0, 3.0]
        ];
        final LU xlu = lu(a);

        expect(xlu.p.matmul(xlu.l).matmul(xlu.u), a);
        expect(xlu.isSingular, false);
      });

      test('3', () {
        final a = [
          [1.0, 2.0, 4.0],
          [3.0, 8.0, 14.0],
          [2.0, 6.0, 13.0]
        ];
        final LU xlu = lu(a);

        expect(xlu.p.matmul(xlu.l).matmul(xlu.u), a);
        expect(xlu.isSingular, false);
      });

      test('4', () {
        final a = [
          [2.0, 5.0, 3.0, 5.0],
          [4.0, 6.0, 6.0, 3.0],
          [11.0, 3.0, 2.0, -2.0],
          [4.0, -7.0, 9.0, 3.0]
        ];
        final LU xlu = lu(a);

        expect(xlu.p.matmul(xlu.l).matmul(xlu.u).isClose(a), isTrue);
        expect(xlu.isSingular, false);
      });
    });

    group(('Determinant'), () {
      test('NotSquare', () {
        final a = [
          [6.0, 1.0],
          [4.0, -2.0],
          [2.0, 8.0]
        ];
        final LU xlu = lu(a);

        expect(xlu.p.matmul(xlu.l).matmul(xlu.u), a);
        expect(() => xlu.determinant, throwsUnsupportedError);
      });

      test('0', () {
        final a = [
          [3.0, 8.0],
          [4.0, 6.0]
        ];
        final LU xlu = lu(a);

        expect(xlu.p.matmul(xlu.l).matmul(xlu.u).isClose(a), isTrue);
        expect(xlu.determinant, -14);
      });

      test('1', () {
        final a = [
          [6.0, 1.0, 1.0],
          [4.0, -2.0, 5.0],
          [2.0, 8.0, 7.0]
        ];
        final LU xlu = lu(a);

        expect(xlu.p.matmul(xlu.l).matmul(xlu.u).isClose(a), isTrue);
        expect(xlu.determinant, -306);
      });

      test('2', () {
        final a = [
          [2.0, 5.0, 3.0, 5.0],
          [4.0, 6.0, 6.0, 3.0],
          [11.0, 3.0, 2.0, -2.0],
          [4.0, -7.0, 9.0, 3.0]
        ];
        final LU xlu = lu(a);

        expect(xlu.p.matmul(xlu.l).matmul(xlu.u).isClose(a), isTrue);
        expect(xlu.determinant.round(), 2960);
      });
    });

    group('Solve', () {
      test('SizeMismatch', () {
        final a = [
          [2.0, 1.0],
          [0.0, -1.0],
          [-2.0, 3.0],
          [-1.0, 0.0]
        ];
        final b = [
          [0.0, 1.0, 11.0],
          [2.0, -1.0, -5.0],
          [-8.0, 3.0, 9.0]
        ];
        final LU xlu = lu(a);

        expect(() => xlu.solve(b), throwsArgumentError);
      });

      test('Singular', () {
        final a = [
          [0.0, 0.0],
          [0.0, 1.0]
        ];
        final b = [
          [0.0, 1.0, 11.0],
          [2.0, -1.0, -5.0]
        ];
        final LU xlu = lu(a);

        expect(() => xlu.solve(b), throwsUnsupportedError);
      });

      test('Identity', () {
        final a = [
          [6.0, 1.0, 1.0],
          [4.0, -2.0, 5.0],
          [2.0, 8.0, 7.0]
        ];
        final LU xlu = lu(a);
        Double2D solution = xlu.solve(a);

        expect(solution, [1.0, 1.0, 1.0].diagonal(0));
      });

      test('0', () {
        final a = [
          [3.0, 2.0],
          [-6.0, 6.0]
        ];
        final b = [
          [7.0],
          [6.0]
        ];
        final LU xlu = lu(a);
        Double2D solution = xlu.solve(b);
        Double2D product = a.matmul(solution);

        expect(product, b);
      });

      test('1', () {
        var a = [
          [6.0, 3.0, 0.0],
          [2.0, 5.0, 1.0],
          [9.0, 8.0, 6.0]
        ];
        var b = [
          [60.0, 45.0],
          [49.0, 43.0],
          [141.0, 92.0]
        ];
        final LU xlu = lu(a);
        Double2D solution = xlu.solve(b);
        Double2D product = a.matmul(solution);

        expect(product, b);
      });

      test('2', () {
        final a = [
          [1.0, 4.0],
          [2.0, 5.0],
        ];
        final b = a.matmul([5, 2].toCol());
        final LU xlu = lu(a);

        expect(
            xlu.solve(b),
            [
              [5],
              [2]
            ]);
      });
    });
  });
}
