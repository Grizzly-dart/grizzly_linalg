import 'package:grizzly/grizzly.dart';
import 'package:grizzly_linalg/grizzly_linalg.dart';

main() {
  Double2D a = [
    [1, 1],
    [1, 2],
    [1, 3],
    [1, 4],
    [1, 5],
  ].toDouble();

  final asvd = svd(a);

  Double2D u = asvd.u;
  Double1D s = asvd.s;
  Double2D v = asvd.v;

  print('U');
  print(u);
  print('S');
  print(s);
  print('V');
  print(v);

  print('A = U * S * V^T');
  print(u.matmulDiag(s).matmul(v.transpose));

  Double2D b = (a.cols[0].plus(a.cols[1] * 5)).toCol();
  print('B');
  print(b);

  print('x = A^-1 * b');
  print(asvd.solve(b));
}
