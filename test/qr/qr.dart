import 'package:grizzly_array/grizzly_array.dart';
import 'package:grizzly_linalg/grizzly_linalg.dart';
import 'package:test/test.dart';

void main() {
  group('QR', () {
    group('iSFullRank', () {
      test('0', () {
        final a = new Double2D([
          [1.0, 0.0, 2.0],
          [2.0, 1.0, 0.0],
          [3.0, 2.0, 1.0]
        ]);
        QR aqr = qr(a);

        expect(aqr.isFullRank, isTrue);
      });

      test('1', () {
        final a = new Double2D([
          [2.0, 5.0, 3.0],
          [4.0, 6.0, 6.0],
          [11.0, 3.0, 2.0],
          [4.0, -7.0, 9.0]
        ]);
        QR aqr = qr(a);

        expect(aqr.isFullRank, isTrue);
      });

      test('2', () {
        final a = new Double2D([
          [1.0, 0.0, 1.0],
          [2.0, 1.0, 3.0],
          [3.0, 2.0, 5.0]
        ]);
        QR aqr = qr(a);

        expect(aqr.isFullRank, isFalse);
      });

      test('3', () {
        final a = new Double2D([
          [1.0, 2.0],
          [2.0, 4.0],
          [3.0, 6.0]
        ]);
        QR aqr = qr(a);

        expect(aqr.isFullRank, isFalse);
      });

      test('4', () {
        final a = new Double2D([
          [1.0, 0.0, 2.0],
          [2.0, 1.0, 0.0]
        ]);
        QR aqr = qr(a);

        expect(aqr.isFullRank, isFalse);
      });
    });

    group('Factorization', () {
      test('0', () {
        final a = new Double2D([
          [1.0, 1.0],
          [1.0, 2.0],
          [1.0, 3.0],
          [1.0, 4.0],
          [1.0, 5.0],
        ]);

        final aqr = qr(a);
        expect(
            aqr.q,
            doubles2([
              [
                -0.44721359549995787,
                -0.632455532033676,
              ],
              [
                -0.4472135954999579,
                -0.3162277660168379,
              ],
              [
                -0.4472135954999579,
                8.326672684688674e-17,
              ],
              [
                -0.4472135954999579,
                0.3162277660168379,
              ],
              [
                -0.4472135954999579,
                0.6324555320336758,
              ],
            ]));
        expect(
            aqr.r,
            doubles2([
              [
                -2.23606797749979,
                -6.708203932499368,
              ],
              [
                0.0,
                3.16227766016838,
              ],
            ]));

        expect(aqr.q.matmul(aqr.r).isClose(a), true);
      });

      test('1', () {
        final a = new Double2D([
          [1.0],
          [2.0],
          [3.0],
          [4.0],
          [5.0],
        ]);

        final aqr = qr(a);
        expect(
            aqr.q,
            doubles2([
              [
                -0.13483997249264834,
              ],
              [
                -0.26967994498529685,
              ],
              [
                -0.40451991747794525,
              ],
              [
                -0.5393598899705937,
              ],
              [
                -0.674199862463242,
              ],
            ]));
        expect(
            aqr.r,
            doubles2([
              [
                -7.416198487095663,
              ],
            ]));

        expect(aqr.q.matmul(aqr.r).isClose(a), true);
      });

      test('2', () {
        final a = new Double2D([
          [2.0, 5.0, 3.0],
          [4.0, 6.0, 6.0],
          [11.0, 3.0, 2.0],
          [4.0, -7.0, 9.0]
        ]);
        final aqr = qr(a);

        expect(
            aqr.q,
            doubles2([
              [
                -0.1596173768935245,
                -0.4307106780299486,
                0.35286080004698334,
              ],
              [
                -0.31923475378704885,
                -0.47883817953541685,
                0.5868695178973314,
              ],
              [
                -0.8778955729143844,
                -0.02558677295227413,
                -0.4776908104748165,
              ],
              [
                -0.31923475378704885,
                0.764557144169145,
                0.5503498108849223,
              ],
            ]));

        expect(
            aqr.r,
            doubles2([
              [
                -12.529964086141668,
                -3.1125388494237267,
                -7.023164583315076,
              ],
              [
                0.0,
                -10.45524279540308,
                2.66467964031541,
              ],
              [
                0.0,
                0.0,
                8.577566184539606,
              ],
            ]));

        expect(aqr.q.matmul(aqr.r).isClose(a), true);
      });
    });

    group('Solve', () {
      test('InvalidSize', () {
        final a = new Double2D([
          [2.0, 1.0],
          [0.0, -1.0],
          [-2.0, 3.0],
          [-1.0, 0.0]
        ]);
        final aqr = qr(a);
        final b = new Double2D([
          [0.0, 1.0, 11.0],
          [2.0, -1.0, -5.0],
          [-8.0, 3.0, 9.0]
        ]);

        expect(() => aqr.solve(b), throwsArgumentError);
      });

      test('NotFullRank', () {
        final a = new Double2D([
          [1.0, 2.0],
          [2.0, 4.0],
          [3.0, 6.0]
        ]);
        final aqr = qr(a);
        final b = new Double2D([
          [0.0, 1.0],
          [2.0, -1.0],
          [8.0, 11.0]
        ]);

        expect(() => aqr.solve(b), throwsUnsupportedError);
      });

      test('0', () {
        final a = new Double2D([
          [1.0, 1.0],
          [1.0, 2.0],
          [1.0, 3.0],
          [1.0, 4.0],
          [1.0, 5.0],
        ]);
        final aqr = qr(a);
        final b = a.matmul(new Double2D.aCol(<int>[5, 2]));

        expect(
            aqr.solve(b).isClose(doubles2([
                  [5.0],
                  [2.0]
                ])),
            true);
      });

      test('Identity', () {
        final a = new Double2D([
          [6.0, 1.0],
          [4.0, -2.0],
          [2.0, 8.0]
        ]);
        final aqr = qr(a);
        Double2D solution = aqr.solve(a);

        expect(solution.isClose(doubles([1, 1]).diagonal()), isTrue);
      });

      test('NotFullRank', () {
        final a = new Double2D([
          [4.0, 8.0],
          [0.0, 2.0],
          [1.0, 6.0]
        ]);
        final aqr = qr(a);
        final b = new Double2D([
          [92.0, 40.0],
          [18.0, 8.0],
          [59.0, 26.0]
        ]);
        Double2D solution = aqr.solve(b);
        final product = a.matmul(solution);

        expect(b.isClose(product), isTrue);
      });
    });
  });
}
