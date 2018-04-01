import 'package:grizzly_array/grizzly_array.dart';
import 'package:grizzly_linalg/grizzly_linalg.dart';
import 'package:test/test.dart';

void main() {
  group('SVD', () {
    group('Factorization', () {
      test('0', () {
        Double2D x = new Int2D([
          [1, 1, 5],
          [1, 2, 8],
          [1, 3, 4],
          [1, 4, 10],
          [1, 5, 11],
        ]).toDouble();

        final xsvd = svd(x);

        expect(
            xsvd.u,
            doubles2([
              [
                -0.26087531294877375,
                -0.5684718823447511,
                -0.4594088781793647,
              ],
              [
                -0.421608500321257,
                -0.10734577252718541,
                -0.5429678320750794,
              ],
              [
                -0.2511873639245863,
                -0.6813984899936779,
                0.6395034342383077,
              ],
              [
                -0.5538438329124213,
                0.2233756886599124,
                0.013360100252978098,
              ],
              [
                -0.6199614992080036,
                0.38873641925346136,
                0.2915240664170067,
              ],
            ]));

        expect(
            xsvd.s,
            doubles([
              19.52774872534074,
              0.7497923001663147,
              2.02604077613584,
            ]));

        expect(
            xsvd.v,
            doubles2([
              [
                -0.10792214396840416,
                -0.9937472507879407,
                -0.028621886602278317,
              ],
              [
                -0.3673151131726066,
                0.013104900549950463,
                0.9300042307519741,
              ],
              [
                -0.9238140605534527,
                0.1108813019989689,
                -0.36643269285198243,
              ],
            ]));

        expect(xsvd.u.matmulDiag(xsvd.s).matmul(xsvd.v.transpose).isClose(x),
            true);
      });
    });

    group('Solve', () {
      test('0', () {
        Double2D a = new Int2D([
          [1, 1],
          [1, 2],
          [1, 3],
          [1, 4],
          [1, 5],
        ]).toDouble();
        Double2D b = new Double2D.aCol(a.col[0] + (a.col[1] * 5));

        final asvd = svd(a);
        expect(
            asvd.solve(b),
            doubles2([
              [1.0],
              [5.0]
            ]));
      });
    });
  });
}
