#!/usr/bin/env k8

Math.solve = function(a, b) { // Gauss-Jordan elimination, translated from gaussj() in Numerical Recipes in C.
	// on return, a[n][n] is the inverse; b[n][m] is the solution
	var n = a.length, m = (b)? b[0].length : 0;
	if (a[0].length != n || (b && b.length != n)) return -1; // invalid input
	var xc = [], xr = [], ipiv = [];
	var i, ic, ir, j, l, tmp;

	for (j = 0; j < n; ++j) ipiv[j] = 0;
	for (i = 0; i < n; ++i) {
		var big = 0;
		for (j = 0; j < n; ++j) {
			if (ipiv[j] != 1) {
				for (k = 0; k < n; ++k) {
					if (ipiv[k] == 0) {
						if (Math.abs(a[j][k]) >= big) {
							big = Math.abs(a[j][k]);
							ir = j; ic = k;
						}
					} else if (ipiv[k] > 1) return -2; // singular matrix
				}
			}
		}
		++ipiv[ic];
		if (ir != ic) {
			for (l = 0; l < n; ++l) tmp = a[ir][l], a[ir][l] = a[ic][l], a[ic][l] = tmp;
			if (b) for (l = 0; l < m; ++l) tmp = b[ir][l], b[ir][l] = b[ic][l], b[ic][l] = tmp;
		}
		xr[i] = ir; xc[i] = ic;
		if (a[ic][ic] == 0) return -3; // singular matrix
		var pivinv = 1. / a[ic][ic];
		a[ic][ic] = 1.;
		for (l = 0; l < n; ++l) a[ic][l] *= pivinv;
		if (b) for (l = 0; l < m; ++l) b[ic][l] *= pivinv;
		for (var ll = 0; ll < n; ++ll) {
			if (ll != ic) {
				var dum = a[ll][ic];
				a[ll][ic] = 0;
				for (l = 0; l < n; ++l) a[ll][l] -= a[ic][l] * dum;
				if (b) for (l = 0; l < m; ++l) b[ll][l] -= b[ic][l] * dum;
			}
		}
	}
	for (l = n - 1; l >= 0; --l)
		if (xr[l] != xc[l])
			for (var k = 0; k < n; ++k)
				tmp = a[k][xr[l]], a[k][xr[l]] = a[k][xc[l]], a[k][xc[l]] = tmp;
	return 0;
}

var buf = new Bytes();
var file = arguments.length > 0? new File(arguments[0]) : new File();

var a = []
while (file.readline(buf) >= 0) {
	var t = buf.toString().split("\t");
	if (t[0] != 'CT') continue;
	a[parseInt(t[1])] = [parseInt(t[2]), parseInt(t[3])];
}

file.close();
buf.destroy();

var sum_sr = 0, sum_asm = 0;
for (var i = 0; i < a.length; ++i)
	sum_sr += a[i][0], sum_asm += a[i][1];
var max_q = 0, max_cnt = 0;
for (var i = 0; i < a.length - 1; ++i)
	if (max_cnt < a[i][1]) max_cnt = a[i][1], max_q = i;
var r = a[max_q][1] / a[max_q][0];

var f = 0.00008;
var b = [];
b[0] = [a[0][0], a[0][1]];
b[1] = [a[1][0], a[1][1]];
for (var i = 2; i < a.length; ++i) {
	if (i <= max_q) {
		var y = (r * a[i][0] - a[i][1]) / (r - f);
		var z = a[i][0] - y;
		if (z < 0) z = 0;
		b[i] = [z, r * z];
	} else b[i] = [a[i][0], a[i][1]];
	//print(i, a[i][0], a[i][1], b[i][0].toFixed(0), b[i][1].toFixed(0));
}

var min_cnt = max_cnt, min_q = max_q;
for (var i = max_q; i >= 2; --i)
	if (min_cnt > a[i][1]) min_cnt = a[i][1], min_q = i;

if (max_cnt - min_cnt + 1 < 5) throw("ERROR: not enough points");

var N = 2;
var x = [], y = [], cap = max_q;
if (max_q > min_q + 8) cap = min_q + 8;
for (var i = min_q; i < cap; ++i)
	x[i - min_q] = i, y[i - min_q] = b[i+1][1] / b[i][1];
var xk = [];
for (var i = 0; i < x.length; ++i) {
	var t = 1;
	xk[i] = [];
	for (var k = 0; k <= N * N; ++k) {
		xk[i][k] = t;
		t *= x[i];
	}
	//print(x[i], y[i]);
}

var A = [], B = [];
for (var i = 0; i <= N; ++i) A[i] = [];
for (var i = 0; i <= N; ++i) {
	for (var j = 0; j <= N; ++j) {
		var sum = 0;
		for (var l = 0; l < x.length; ++l)
			sum += xk[l][i + j];
		A[i][j] = A[j][i] = sum;
	}
	var sum = 0;
	for (var l = 0; l < x.length; ++l)
		sum += xk[l][i] * y[l];
	B[i] = [sum];
}

Math.solve(A, B);

function poly_cal(B, x)
{
	var t = 1, s = 0;
	for (var i = 0; i < B.length; ++i) {
		s += t * B[i][0];
		t *= x;
	}
	return s;
}

//print(B.join("\t"));
//print(poly_cal(B, 1), poly_cal(B, 0));
for (var i = min_q - 1; i >= 0; --i) {
	var r = poly_cal(B, i);
	if (r < 1) r = 1;
	b[i][1] = b[i+1][1] / r;
}

//for (var i = 0; i < b.length; ++i) print(i, b[i].join("\t"));

var sum_asm2 = 0;
for (var i = 0; i < b.length; ++i)
	sum_asm2 += b[i][1];

print(sum_asm - sum_asm2, -4.343 * Math.log(Math.log(sum_asm / sum_asm2) / 31));
