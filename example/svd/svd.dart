import 'package:grizzly_array/grizzly_array.dart';
import 'package:grizzly_linalg/grizzly_linalg.dart';

main() {
  Double2D a = new Int2D([
    [1, 1],
    [1, 2],
    [1, 3],
    [1, 4],
    [1, 5],
  ]).toDouble();

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

  Double2D b = new Double2D.aCol(a.col[0] + (a.col[1] * 5));
  print('B');
  print(b);

  print('x = A^-1 * b');
  print(asvd.solve(b));
}
