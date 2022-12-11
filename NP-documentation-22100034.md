# **Numerical programming document**



## Errors and Taylor series



###  Factorial

```c
double factorial(int N)
{
	int y = 1;
	for (int k = 2; k <= N; k++)
		y = y * k;

	return y;
}
```

##  

### sinetaylor ###

Taylor series approximation for sin(x) using pre-defined functions (input unit: [rad])

```c++
double sinTaylor(double _x)
{
	int N_max = 10;
	double S_N = 0;
    for (int k = 0; k < N_max; k++)
	S_N = S_N + pow(-1, k) * pow(_x, 2 * k + 1) / factorial(2 * k + 1);

return S_N;
}
```


### non linear solve 

#### bisection

```c++
double bisection(float _a0, float _b0, float _tol);

Define Bisection function, assuming (func(a) * func(b) <0 )

1. Initialization of K, Nmax, range a,b , xn, epsilon
2. update range a,b 

if (func(xn) * func(a) < 0)
			b = xn;

else a = xn;

3. check the tolerance

 ep = fabs(func(xn))
```



#### NewtonRaphson

```c++
double newtonRaphson(double func(double x), double dfunc(double x), double x0, double tol)
{
	double xn = x0; // initial value
	double ep = 1000; // tolerance error
	int Nmax = 1000; // repeat
	int k = 0;
	double h = 0; 

​    h = -func(xn) / dfunc(xn); // h definition


xn = xn + h; // update xn
ep = fabs(func(xn));


}
```





## Differentiation

### gradient1D

:  Return the dy/dx results for the input data. 

```c++
void gradient1D(double x[], double y[], double dydx[], int m) {
	// Check if length x and y are equal	
double h = 0.02;
// For starting point	
dydx[0] = (-3 * y[0] + 4 * y[1] - y[2]) / (2 * h);// Three-Point FWD  O(h^2). Need to use first 2 points

// For mid points	
for (int i = 1; i < m - 1; i++) {
	dydx[i] = (y[i + 1] - y[i - 1]) / (2 * h);
}
// Two-Point Central  O(h^2)
// For end point	
dydx[m - 1] = (y[m - 3] - 4 * y[m - 2] + 3 * y[m - 1]) / (2 * h); // Three-Point BWD  O(h^2). Need to use last 2 points
}

// Return the d2y/d2x results for the input data. 
```




```c++
void acceleration(double x[], double y[], double dy2dx2[], int m) {
double h = x[1] - x[0]; // define h

dy2dx2[0] = (2 * y[0] - 5 * y[1] + 4 * y[2] - y[3]) / (pow(h, 2));// for starting point, 4-point FWD

for (int i = 1; i < m - 1; i++) {
	dy2dx2[i] = (y[i - 1] - 2 * y[i] + y[i + 1]) / (pow(h, 2)); // for mid point, 3-point central
}

dy2dx2[m - 1] = (-y[m - 4] + 4 * y[m - 3] - 5 * y[m - 2] + 2 * y[m - 1]) / (pow(h, 2)); // for end point, 4-point BWD
```
}

## integratoin ##

### trapzodial ###

```c++
double trapz(double x[], double y[], int m) {
int N = m - 1; // N is interval number
double I = 0;

for (int i = 0; i < N; i++) {
	I += (y[i] + y[i + 1]) * (x[i + 1] - x[i]);
}
return I * 0.5;
}
```


### simpson ###

```c++
double simpson13(double x[], double y[], int m) {
	double h = 0.75;
	double I = y[0] + 4 * y[m - 2] + 4 * y[m - 1];
	for (int i = 1; i < m - 3; i += 2) {
		I += 4 * y[i] + 2 * y[i + 1];
	}
	return I * h / 3;
}
```

### integral ###

```c++
double integral(double func(const double x), double a, double b, int n) {
double h = (b - a) / n; // a is x[0], b is x[N]
double I = func(a) + 4 * func(b - h) + func(b);
for (int i = 1; i < n - 2; i += 2) {
	// x0 = a;
	// xi = a+i*h;
	double xi = a + i * h;
	I += 4 * func(xi) + 2 * func(xi + h);
}
return (I * h / 3);
}
```
### steffensen ###



```c++
double steffensen(double func(double x), double x0, double tol) {
	double sol = 0;
	//double tol = 1 * pow(10, -5);
	double xn = x0;
	double ep = 1000;
	int Nmax = 1000;
	int k = 0;
	double h = 0;
	double g = 0;
do {
	if (dfunc(xn) == 0)
	{
		printf("[ERROR] dF == 0 !!\n");
		break;
	}
	else
	{
		g = func(xn + func(xn));

		h = -func(xn) * func(xn) / (g - func(xn)); // h definition

		xn = xn + h; // update xn

		ep = fabs(func(xn));

		k++;

		printf("k:%d \t", k);
		printf("X(k): %f \t", xn);
		printf("Tol: %.10f\n", ep);

	}
} while (k < Nmax && ep > tol);

return xn;
}
```


## linear solve ## 

### gauss elimination ###

void gaussElim(Matrix _A, Matrix _B, Matrix* _U, Matrix* _B_out);

**Parameters**

- **A**: Matrix **A** in structure Matrix form. Should be (nxn) square.

- **B**: vector **b** in structure Matrix form. Should be (nx1)

- **U**: Matrix **U** in structure Matrix form. Should be (nxn) square.

- **B_out**: vector **B_out** in structure Matrix form. Should be (nx1)

  

**Example code**

``` c++
Matrix prob1_matA = txt2Mat(path, "prob1_matA");
Matrix prob1_vecb = txt2Mat(path, "prob1_vecb");
Matrix prob1_matU = createMat(prob1_matA.rows, prob1_matA.cols);
Matrix prob1_vecd = createMat(prob1_matA.rows, 1);
Matrix prob1_vecx = createMat(prob1_vecb.rows, 1);

gaussElim(prob1_matA, prob1_vecb, prob1_matU, prob1_vecd);
backSub(prob1_matU, prob1_vecd, prob1_vecx);


void gaussElim(Matrix _A, Matrix _b, Matrix _U, Matrix _d) {
	double s = 0;
	int m = _A.rows; // number of rows
	int n = _A.cols;
	int i, j, k;

	for (int i = 0; i < _A.rows; i++) {
		for (int j = 0; j < _A.cols; j++) {
			_U.at[i][j] = _A.at[i][j];
		}
	}
	for (int i = 0; i < _b.rows; i++) {
		for (int j = 0; j < _b.cols; j++) {
			_d.at[i][j] = _b.at[i][j];
		}
	}

	for (k = 0; k < m - 1; k++) {
		for (i = k + 1; i < m; i++) {
			s = _U.at[i][k] / _U.at[k][k];
			_U.at[i][k] = 0.0;
			for (j = k + 1; j < n; j++) {
				_U.at[i][j] = _U.at[i][j] - s * _U.at[k][j];
			}
			_d.at[i][0] = _d.at[i][0] - s * _d.at[k][0];
		}
	}


} 
```



### LU decomposition ###

void LUdecomp(Matrix A, Matrix L, Matrix U);

**Parameters**

- **A**: Matrix **A** in structure Matrix form. Should be (nxn) square.

- **L**: Matrix **L** in structure Matrix form. Should be (nxn) square. Lower triangular matrix.
- **U**: Matrix **A** in structure Matrix form. Should be (nxn) square. Upper triangular matrix.

``` c++
Matrix prob1_matA = txt2Mat(path, "prob1_matA");
Matrix prob1_vecb = txt2Mat(path, "prob1_vecb");
Matrix prob1_matL = zeros(prob1_matA.rows, prob1_matA.cols);
Matrix prob1_matU = zeros(prob1_matA.rows, prob1_matA.cols);

LUdecomp(prob1_matA, prob1_matL, prob1_matU);
solveLU(prob1_matL, prob1_matU, prob1_vecb, prob1_vecx);
invMat(prob1_matA, Ainv);
x_check = multiplyMat(Ainv, prob1_vecb);

void LUdecomp(Matrix A, Matrix L, Matrix U) {
	double m = 0;
	// A -> U
	for (int i = 0; i < A.rows; i++) {
		for (int j = 0; j < A.cols; j++) {
			U.at[i][j] = A.at[i][j];
		}
	}
	// LU
	for (int k = 0; k <= U.rows-2; k++) {
		for (int i = k+1; i<= U.rows-1; i++) {

			m = U.at[i][k] / U.at[k][k];
			L.at[i][k] = m;

			for (int j = 0; j < U.rows; j++) {
				U.at[i][j] = U.at[i][j] - m * U.at[k][j];

			}

		}
	}

	for (int i = 0; i < L.rows; i++) {
		L.at[i][i] = 1;
	}
}
```



### LU decomposition with Pivoting ###

void LUdecomp(Matrix A, Matrix L, Matrix U, Matrix P);

**Parameters**

- **A**: Matrix **A** in structure Matrix form. Should be (nxn) square.

- **L**: Matrix **L** in structure Matrix form. Should be (nxn) square. Lower triangular matrix.
- **U**: Matrix **A** in structure Matrix form. Should be (nxn) square. Upper triangular matrix.
- **P**: Matrix **P** in structure Matrix form. Should be (nxn) square. Pivoting matirx.

``` c++
Matrix prob1P_vecb = txt2Mat(path, "prob1_vecb");
Matrix prob1_vecx = createMat(prob1_vecb.rows, 1);
Matrix prob1P_vecx = createMat(prob1_vecb.rows, 1);
Matrix Ainv = createMat(prob1_matA.rows, prob1_matA.cols);
Matrix prob1P_matP = eye(prob1_matA.rows, prob1_matA.cols);
Matrix prob1P_matL = zeros(prob1_matA.rows, prob1_matA.cols);
Matrix prob1P_matU = zeros(prob1_matA.rows, prob1_matA.cols);
Matrix x_check = createMat(prob1_matA.rows, 1);

LUdecomp(prob1_matA, prob1P_matL, prob1P_matU, prob1P_matP);
solveLU(prob1P_matL, prob1P_matU, prob1P_matP, prob1P_vecb, prob1P_vecx);
	
void LUdecomp(Matrix A, Matrix L, Matrix U, Matrix P) {
	Matrix s = zeros(A.rows, 1);
	double sp = 0;
	double rest_U = 0;
	double rest_P = 0;
	double rest_L = 0;
	double m = 0.0;

	// defind U = A
	for (int i = 0; i < A.rows; i++) {
		for (int j = 0; j < A.cols; j++) {
			U.at[i][j] = A.at[i][j];
		}
	}

	for (int k = 0; k < U.rows-1; k++) {
		// find s which is the biggest number of each rows
		double max = 0.0;
		for (int i = k; i < U.rows; i++) {
			max = U.at[i][i];
			for (int j = 0; j < A.cols; j++) {
				if (fabs(U.at[i][j]) > fabs(max)) {
					max = U.at[i][j];

				}
			} s.at[i][0] = max;
		}
;		// find spmax
		int spmax_row = 0;
		double spmax = 0;

		for (int i = k; i < U.rows; i++) {

				sp = fabs(U.at[i][k]) / fabs(s.at[i][0]);
				if (sp > spmax) {
					spmax = sp;
					spmax_row = i;
				}
		}

		//row interchange for U,L,P
		
			for (int j = 0; j < U.cols; j++) {

				// A row interchange
				rest_U = U.at[k][j];
				U.at[k][j] = U.at[spmax_row][j];
				U.at[spmax_row][j] = rest_U;
				//L row interchage
				rest_L = L.at[k][j];
				L.at[k][j] = L.at[spmax_row][j];
				L.at[spmax_row][j] = rest_L;
				// P row interchange
				rest_P = P.at[k][j];
				P.at[k][j] = P.at[spmax_row][j];
				P.at[spmax_row][j] = rest_P;

			}

		// find m to make matrix L
		for (int i = k + 1; i < U.rows ; i++) {
		
			m = U.at[i][k] / U.at[k][k];
			L.at[i][k] = m;

			for (int j = k; j < U.cols; j++) {
				U.at[i][j] = U.at[i][j] - m * U.at[k][j];

			}
		}
		
	} // L = L + I
	for (int i = 0; i < L.rows; i++) {
		L.at[i][i] = 1;
	} 
}

```

```c++
void solveLU(Matrix L, Matrix U, Matrix P, Matrix b, Matrix x) {
	// update b = Pb

	Matrix bP = zeros(b.rows, b.cols);


	for (int i = 0; i < P.rows; i++) {
		for (int j = 0; j < b.cols; j++) {
			for (int k = 0; k < P.rows; k++) {
				bP.at[i][j] =bP.at[i][j]+ P.at[i][k] * b.at[k][j];
			}
		} 
	}
	Matrix y = zeros(L.rows,1);
	fwdsub(L, bP, y);
	backsub(U, y, x);

}
void solveLU(Matrix L, Matrix U, Matrix b, Matrix x) {

	Matrix y = zeros(b.rows, 1);
	fwdsub(L, b, y);
	backsub(U, y, x);
}
```



### Norm ###

- **vnorm** is the scale of matrix

```c++
double vnorm(Matrix x) {
	double s = 0;
	for (int i = 0; i < x.rows; i++) { 
		for (int j = 0; j < x.cols; j++) {
			s += x.at[i][j] * x.at[i][j];
		}
	}return sqrt(s);
}
```



### Inverse Matrix ### 

void invMat(Matrix A, Matrix Ainv);

**Parameters**

- **A**: Matrix **A** in structure Matrix form. Should be (nxn) square.

- **Ainv**: Inverse Matrix of  **A** in structure Matrix form. Should be (nxn) square.

```c++
void invMat(Matrix A, Matrix Ainv) {
	Matrix L = zeros(A.rows, A.cols);
	Matrix U = zeros(A.rows, A.cols);
	Matrix L_inverse = createMat(L.rows, 1);
	Matrix a_inverse = createMat(A.rows, 1);
	Matrix A_inverse = createMat(A.rows, A.cols);
	Matrix I = eye(A.rows, A.cols);
	Matrix c = zeros(A.rows, 1);

	LUdecomp(A, L, U);


	for (int j = 0; j < A.cols; j++) {
		for (int i = 0; i < A.rows; i++) {
			c.at[i][0] = I.at[i][j];
		}
		fwdsub(L, c, L_inverse); // to get L_inverse
		backsub(U, L_inverse, a_inverse); // to get A_inverse


		// update A_inverse 
		for (int k = 0; k < A.rows; k++) {
			A_inverse.at[k][j] = a_inverse.at[k][0];
		}
	} 

	for (int i = 0; i < A.rows; i++) {
		for (int j = 0; j < A.cols; j++) {
			Ainv.at[i][j] = A_inverse.at[i][j];
		}

	}
}
```



### Eigenvalue  ###

: find eigenvalue through QR decomposition from matrix A

: by updated U = Q * R, we can find eigenvalue. 

Matrix eig(Matrix A);

**Parameters**

- **A**: Matrix **A** in structure Matrix form. Should be (nxn) square.

```c++
Matrix matA1 = txt2Mat(path, "prob_matA");
eig(matA1);	// print matirx in function. don't need to print matrix in main loop

Matrix eig(Matrix A) {
	Matrix U = createMat(A.rows, A.cols); // U = A
	Matrix Q = eye(A.rows, A.cols);
	Matrix R = copyMat(A);
	Matrix lamda = createMat(A.rows, 1);

	if (A.rows != A.cols) {
		printf("\n****************************************************");
		printf("\n  ERROR!!: A should be n by b matrix");
		printf("\n****************************************************\n");
		return createMat(0, 0);
	}
	else {
        
		// U = A, initialize U
		Matrix U = copyMat(A);
		//printMat(U, "UUUUUUUUUUU");
		
		// step1. A = QR using HouseHold Matrix
		for (int i = 0; i <500; i++) {
			QRdecomp(U,Q,R);
			
		// step2. A = RQ , make into similar matrix
			U = multiplyMat(R, Q);
			//printMat(Q, "QQQQQQQQQQQQq");
		}
		
		lamda = diag(U); // find eigenvalue 
		printMat(U, "final matrix U:\n");
		printMat(lamda, "Eigenvaluse of A::\n");
	}
	
	return lamda;
}
```

### Eigenvector ###

: with eigenvalue, we can find eigenvector easily by simple method

: B = A - lamda * I // lamda is eigenvalue and  I is identity matrix

<img src="C:\Users\eunji\AppData\Roaming\Typora\typora-user-images\image-20221114192455278.png" alt="image-20221114192455278" style="zoom: 67%;" />

Matrix eigvec(Matrix A);

**Parameters**

- **A**: Matrix **A** in structure Matrix form. Should be (2 by 2 or 3 by 3) square.

```c++
Matrix matA2 = txt2Mat(path, "prob_matA_eigvec");
eigvec(matA2); // print matirx in function. don't need to print matrix in main loop

Matrix eigvec(Matrix A) {
	Matrix v1 = createMat(A.rows - 1, 1);
	Matrix v2 = createMat(A.rows - 1, 1);
	Matrix v3 = createMat(A.rows - 1, 1);
	Matrix V1 = ones(A.rows, 1);
	Matrix V2 = ones(A.rows, 1);
	Matrix V3 = ones(A.rows, 1);
	Matrix eigvec = zeros(A.rows, A.cols);
	//printMat(A, "AAAAAAAAaa");


	if (A.rows == 3 && A.cols == 3) {
		Matrix lamda = eig(A);
		//printMat(lamda, "lamdaaaaaaaaaaaaa");
		double lamda1 = lamda.at[0][0];
		double lamda2 = lamda.at[1][0];
		double lamda3 = lamda.at[2][0];
		//printf("%f\n", lamda1);
		//printf("%f\n", lamda2);
		//printf("%f\n", lamda3);

		Matrix I = eye(A.rows, A.cols);
		Matrix B1 = addMat(A, mulConstMat(-lamda1, I));
		Matrix B2 = addMat(A, mulConstMat(-lamda2, I));
		Matrix B3 = addMat(A, mulConstMat(-lamda3, I));
		//printMat(B1, "B1\n");		
		//printMat(B2, "B2\n");
		//printMat(B3, "B3\n");


		Matrix b1 = createMat(A.rows - 1, A.cols - 1);
		Matrix b1_ = createMat(A.rows - 1, 1);
		Matrix b1inv = createMat(A.rows - 1, A.cols-1);
		

		Matrix b2 = createMat(A.rows - 1, A.cols - 1);
		Matrix b2_ = createMat(A.rows - 1, 1);
		Matrix b2inv = createMat(A.rows - 1, A.cols-1);
		

		Matrix b3 = createMat(A.rows - 1, A.cols - 1);
		Matrix b3_ = createMat(A.rows - 1, 1);
		Matrix b3inv = createMat(A.rows - 1, A.cols-1);
		
			b1.at[0][0] = B1.at[1][1];
			b1.at[0][1] = B1.at[1][2];
			b1.at[1][0] = B1.at[2][1];
			b1.at[1][1] = B1.at[2][2];
			b1_.at[0][0] = -1.0 * B1.at[1][0];
			b1_.at[1][0] = -1.0 * B1.at[2][0];
			//printMat(b1, "b1\n");
			//printMat(b1_, "b1_");
			invMat(b1, b1inv);
			//printMat(b1, "b111111111111");
			//printMat(b1inv,"b1inv");
			v1 = multiplyMat(b1inv, b1_);

			V1.at[1][0] = v1.at[0][0];
			V1.at[2][0] = v1.at[1][0];

			V1 = mulConstMat(1 / vnorm(V1), V1);
			//printMat(V1, "V1111111111111111\n");
			

			b2.at[0][0] = B2.at[0][0];
			b2.at[0][1] = B2.at[0][2];
			b2.at[1][0] = B2.at[2][0];
			b2.at[1][1] = B2.at[2][2];
			b2_.at[0][0] = -1.0 * B2.at[0][1];
			b2_.at[1][0] = -1.0 * B2.at[2][1];
			//printMat(b2, "b2\n");
			//printMat(b2_, "b2_");
			invMat(b2, b2inv);

			//printMat(b2, "b22222222222222");
			//printMat(b2inv,"b2inv");
			v2 = multiplyMat(b2inv, b2_);

			V2.at[0][0] = v2.at[0][0];
			V2.at[2][0] = v2.at[1][0];

			V2 = mulConstMat(1 / vnorm(V2), V2);
			//printMat(V2, "V2222222222222\n");
			

			b3.at[0][0] = B3.at[0][0];
			b3.at[0][1] = B3.at[0][1];
			b3.at[1][0] = B3.at[1][0];
			b3.at[1][1] = B3.at[1][1];
			b3_.at[0][0] = -1.0 * B3.at[0][2];
			b3_.at[1][0] = -1.0 * B3.at[1][2];
			//printMat(b3, "b3\n");
			//printMat(b3_, "b3_");
			invMat(b3, b3inv);
			//printMat(b3, "b333333333333");
			//printMat(b3inv,"b3inv");

			v3 = multiplyMat(b3inv, b3_);

			
			V3.at[0][0] = v3.at[0][0];
			V3.at[1][0] = v3.at[1][0];
			//printMat(V3,"V33333333");
			
			V3 = mulConstMat(1 / vnorm(V3), V3);
			//printMat(V3, "V333333333333");
			
			
			for (int i = 0; i < A.rows; i++) {
				
				eigvec.at[i][0] = V1.at[i][0];
				eigvec.at[i][1] = V2.at[i][0];
				eigvec.at[i][2] = V3.at[i][0];
				
			}printMat(eigvec, "Eigenvector of A : \n");

	}
	else if(A.rows == 2 && A.cols == 2) {
	Matrix lamda2 = eig(A);
	//printMat(lamda, "lamdaaaaaaaaaaaaa");
	double lamda2_1 = lamda2.at[0][0];
	double lamda2_2 = lamda2.at[1][0];

	//printf("%f\n", lamda1);
	//printf("%f\n", lamda2);

	Matrix I2 = eye(A.rows, A.cols);
	Matrix C1 = addMat(A, mulConstMat(-lamda2_1, I2));
	Matrix C2 = addMat(A, mulConstMat(-lamda2_2, I2));
	//printMat(C1, "C1\n");		
	//printMat(C2, "C2\n");


	double b2_1 = C1.at[1][1];
	double b2_1_ = -1.0 * C1.at[1][0];
	double b2_1inv = 1 / b2_1;

	double b2_2 = C2.at[0][0];
	double b2_2_ = -1.0 * C2.at[0][1];
	double b2_2inv = 1 / b2_2;


	double v1 = b2_1inv * b2_1_;
	double v2 = b2_2inv * b2_2_;

	V1.at[1][0] = v1;
	V1 = mulConstMat(1 / vnorm(V1), V1);
	//printMat(V1, "V1111111111111111\n");
	V2.at[0][0] = v2;
	V2 = mulConstMat(1 / vnorm(V2), V2);
	//printMat(V1, "V1111111111111111\n")Mat(V2, "V2222222222222\n");
	for (int i = 0; i < A.rows; i++) {

		eigvec.at[i][0] = V1.at[i][0];
		eigvec.at[i][1] = V2.at[i][0];

	}printMat(eigvec, "Eigenvector of A : \n");


	}
	
	else {

		printf("\n****************************************************");
		printf("\n  ERROR!!: A should be 3 by 3 matrix or 2 by 2 matrix");
		printf("\n****************************************************\n");
		return createMat(0, 0);
	}
	return eigvec;
}
```

.

### QR decomposition ###

void QRdecomp(Matrix U, Matrix Q, Matrix R);

**Parameters**

- **U**: Matrix **U** in structure Matrix form. Should be (nxn) square.

- **Q**: Matrix **Q** in structure Matrix form. Should be (nxn) square. Lower triangular matrix.
- **R**: Matrix **R** in structure Matrix form. Should be (nxn) square. Upper triangular matrix.

```c++
QRdecomp(U,Q,R);
Matrix U = createMat(A.rows, A.cols); // U = A
Matrix Q = eye(A.rows, A.cols);
Matrix R = copyMat(A);

void QRdecomp(Matrix U, Matrix Q, Matrix R) {

	int n = U.rows;
	Matrix matR = copyMat(U);
	Matrix matQ = eye(n, n);
	Matrix matI = eye(n, n);
	Matrix matH = zeros(n, n);

	Matrix c = zeros(n, 1);
	Matrix v = zeros(n, 1);
	

	for (int j = 0; j < n - 1; j++) {

        // c is n th column of updated R 
		c = extract_cols(R, j);

		for (int i = 0; i < j; i++) {
			c.at[i][0] = 0.0;
		}
		//printMat(c, "ccccccccccccccccccc");
        
		// e is (A.rows,1) column vector with zero and 1 o(or -1)
		Matrix e = zeros(n, 1);

		if (c.at[j][0] >= 0.0) {
			e.at[j][0] = 1.0;
		}
		else {
			e.at[j][0] = -1.0;
		}
		// v = c + norm(c,2)*e
		v = addMat(c, mulConstMat(vnorm(c), e));
		//printf("%f\n", vnorm(c));
		//printMat(v, "vvvvvvvvvvvvvv");
		//printMat(e, "eeeeeeeeeeeee");

		double vtv = 0.0;
		for (int i = 0; i < v.rows; i++) {
			vtv += v.at[i][0] * v.at[i][0];
		}
		//printf("%f\n", vtv);
		

		if (vtv == 0) {
			printf("\n****************************************************");
			printf("\n  ERROR!!: division is zero");
			printf("\n****************************************************\n");
			break;
		}
		else {
			
			Matrix vvt = multiplyMat(v, transpose(v));
			//printMat(vvt, "vvt\n");
			//printMat(transpose(v), "transposev");

			// H = I - 2 * v*vt/vt*v
            matH = addMat(matI, mulConstMat(-2.0 / vtv, vvt));
			
			//printMat(matH, "matH\n");
			//printMat(mulConstMat(-2.0 / vtv, vvt), "HHHHHHHH");
			//printMat(matI, "IIIIIIIii");
            
            // Q (n+1) = Q(n) * H(n)
			matQ = multiplyMat(matQ, matH);
			//printMat(matQ, "QQQQQQQQ");
            
            // R (n+1) = H(n) * R(n) 
			matR = multiplyMat(matH, matR);
			//printMat(matR, "matRRRRRRRRRRRRr");
			
		}
		
	}
	copyVal(matQ, Q);
	copyVal(matR, R);
	//printMat(R, "RRRRRRRRRRRRr");
	//printMat(Q, "QQQQQQQQQQQQ");

}
```





## ODE _ IVP ##

### ODE part1 ##

#### defining function by user _ example  ####

```c++
double myFunc(const double t, const double v);

double myFunc(const double t, const double v) {
	double tau = 1.0;
	double T = 1.0 / tau;
	double f = 10.0;
	double Vm = 1.0;
	double w = 2 * PI * f;

	double dydt = -T * v + 1 * T * Vm * cos(w * t);

	return dydt;
}
```



### Euler's explicit method ###

void odeEU(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, double y0);

: slope of y with initial condition of y(x0) = y0

: y(i+1) = y(i) + slope * h

: x(i+1) = x(i) + h

**Parameters**

- **myfunc**: given function.

- **y[]** : y
- **t0** : initial x
- **tf** : final x
- **h** : x interval
- **y0** : y initial

```c++
int N = 100;
double y[101] = { 0 };
double t0 = 0;
double tf = 0.1;
double h = 0.001;
double y0 = 0;

odeEU(myFunc, y, t0, tf, h, y0);
printVec(y, 101);

void odeEU(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, double y0) {
	
	int N = (tf-t0)/h ; // N =100
	y[0] = y0;
	double t = t0;
	
	for (int i = 0; i < N; i++) {

		t= t0 + h*i;
		//printf("%f\n", t);
		y[i + 1] = y[i] + myfunc(t, y[i]) * h;
		
	} 
	
}
```



### Euler's modified explicit method ###

void odeEM(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, double y0);

: start from the first point of x0, y0, h , estimate of the slope at the end of the interval. 

: using the average of the two slope

slope1 = f(x(i),y(i)) / slope2 = f(x(i+1), yE(i+1))

y(i+1) = y(i) + 1/2 * (slope1 + slope2) * h

repeat until last point

**Parameters**

- **myfunc**: given function.

- **y[]** : y
- **t0** : initial x
- **tf** : final x
- **h** : x interval
- **y0** : y initial

```c++
int N = 100;
double y[101] = { 0 };
double t0 = 0;
double tf = 0.1;
double h = 0.001;
double y0 = 0;

odeEM(myFunc, y, t0, tf, h, y0);
printVec(y, 101);

void odeEM(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, double y0) {
	int N = (tf-t0)/h ; // N =100
	y[0] = y0;
	double t = t0;
	double slope1 = 0;
	double slope2 = 0;
	double yEu = 0;

	for (int i = 0; i < N; i++) {

		t = t0 + h * i;
		//printf("%f\n", t);
		slope1 = myfunc(t, y[i]);
		yEu = y[i] + slope1 * h;
		slope2 = myfunc(t+1, yEu);
	
		y[i + 1] = y[i] + (slope1 + slope2) * h / 2;


	}
}
```



### second order runge- kutta ###

void odeRK2(double myfunc(const double t, const double y), double y[], double t0, double tf,double h, double y0);

y(i+1) = y(i) + C1 * h * f( t , y ) + C2 * h ( f ( t + alpha * h, y + beta * h * f( t,y ) ) )

C1 + C2 = 1 

C2 * alpha = 1/2

C2 * beta = 1/2

// modified euler's method

K1 = f ( t(i), y(i) )

K2 = f(t+alpha, y(i)+ beta * f ( t(i),y(i) ) * h )

f(t,y) = dy / dt

**Parameters**

- **myfunc**: given function.

- **y[]** : y
- **t0** : initial x
- **tf** : final x
- **h** : x interval
- **y0** : y initial

```c++
t0 = 0;
tf = 0.1;
h = 0.001;
y0 = 0;

odeRK2(myFunc, y, t0, tf, h, y0);
printVec(y, 101);

void odeRK2(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, double y0) {
	int N = 100; 
	y[0] = y0;
	double t = 0;
	double C1 = 0.5;
	double C2 = 0.5;

	for (int i = 0; i < N; i++) {
		double K1 = myfunc(t, y[i]);
		double K2 = myfunc(t + h, y[i] + K1 * h);
		t = t0 + h * i;
		//printf("%f\n", t);
		y[i + 1] = y[i] + (C1 * K1 + C2 * K2) * h;

	}
}
```

### Third order runge- kutta ###

void odeRK3(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, double y0);

**Parameters**

- **myfunc**: given function.

- **y[]** : y
- **t0** : initial x
- **tf** : final x
- **h** : x interval
- **y0** : y initial

```c++
t0 = 0;
tf = 0.1;
h = 0.001;
y0 = 0;

odeRK3(myFunc, y, t0, tf, h, y0);
printVec(y, 101);

void odeRK3(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, double y0) {
	int N = (tf-t0)/h; // N =100
	y[0] = y0;
	double t = t0;

	for (int i = 0; i < N; i++) {
		double K1 = myfunc(t, y[i]);
		double K2 = myfunc(t + h/2, y[i] + K1 * h/2);
		double K3 = myfunc(t + h, y[i] - K1 * h + 2 * K2 * h);
		t = t0 + h * i;
		//printf("%f\n", t);
		y[i + 1] = y[i] + (K1 + 4 * K2 + K3) * h / 6;

	}
}
```



### ODE method ###

void ode(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, double y0, int method);

**Parameters**

- **myfunc**: given function.

- **y[]** : y
- **t0** : initial x
- **tf** : final x
- **h** : x interval
- **y0** : y initial

```c++
int method = 0;
	printf("method : ");
	scanf_s("%d", &method);

	ode(myFunc, y, t0, tf, h, y0,method);
	printVec(y, 101);

void ode(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, double y0, int method) {
	if (method == 0) {
		odeEU(myfunc, y, t0, tf, h, y0);

	}
	else if (method == 1) {
		odeRK2(myfunc, y, t0, tf, h, y0);
		
	}
	else if (method == 1) {
		odeRK3(myfunc, y, t0, tf, h, y0);
	
	}
	else return;
}
```



### ODE part 2 ###

#### defining function by user _ example ODE 2

**Parameters**

- **Y[]** : y
- **dYdt[]** : dy/dt
- **t** : x value in for loop 

```c++
void mckfunc(const double t, const double Y[], double dYdt[]);
void mckfunc_free(const double t, const double Y[], double dYdt[]);
void mckfunc_step(const double t, const double Y[], double dYdt[]);

void mckfunc(const double t, const double Y[], double dYdt[])
{
	double m = 1.0;
	double c = 7.0;
	double k = 6.9;
	double f = 5;
	double Fin = 2 * cos(2 * PI * f * t);
	dYdt[0] = Y[1];
	dYdt[1] = (-k * Y[0] - c * Y[1] + Fin) / m; // dYdt = d2y/dt2


}
void mckfunc_free(const double t, const double Y[], double dYdt[])
{
	double m = 1.0;
	double c = 7.0;
	double k = 6.9;
	double f = 5.0;
	double Fin = 0.0;
	dYdt[0] = Y[1];
	dYdt[1] = (-k * Y[0] - c * Y[1] + Fin) / m;



}
void mckfunc_step(const double t, const double Y[], double dYdt[])
{
	double m = 1.0;
	double c = 7.0;
	double k = 6.9;
	double f = 5.0;
	double Fin = 0.0;
	double A = 2.0;

	Fin = A;
	dYdt[0] = Y[1];
	dYdt[1] = (-k * Y[0] - c * Y[1] + Fin) / m;



}

```



### second order runge - kutta with higer order ###

void sys2RK2(void func(const double t, const double Y[], double dYdt[]), double y1[], double y2[], double t0, double tf, double h, double y1_init, double y2_init) ;

**Parameters**

- **myfunc**: given function.

- **y1[]** :  velocity, dy/dt
- **y2[]** : acceleration , d2y/dt2
- **t0** : initial x
- **tf** : final x
- **h** : x interval
- **y1_init** : y (0) initial
- **y2_init** : dydt(0) initial

```c++
double t0 = 0.0;
double tf = 1.0;
double h = 0.01;
double N = (tf - t0) / h;
double y1[101] = { 0 };
double y2[101] = { 0 };
double y1_init = 0.0;
double y2_init = 0.2;

sys2RK2(mckfunc, y1,y2,t0, tf,h,y1_init, y2_init);
sys2RK2(mckfunc_free, y1, y2, t0, tf, h, y1_init, y2_init);
sys2RK2(mckfunc_step, y1, y2, t0, tf, h, y1_init, y2_init);

void sys2RK2(void func(const double t, const double Y[], double dYdt[]), double y1[], double y2[], double t0, double tf, double h, double y1_init, double y2_init) {
	
	double N = (tf-t0)/h;
	double t[101] = { 0.0 };
	double K1[2] = { 0.0 };
	double K2[2] = { 0.0 };
	double Yin_k1[2] = { 0.0 };
	double Yin_k2[2] = { 0.0 };
	double K1_y = 0.0;
	double K1_z = 0.0;
	double K2_y = 0.0;
	double K2_z = 0.0;
	double C1 = 0.5;
	double C2 = 0.5;

	t[0] = t0;
	y1[0] = y1_init;
	y2[0] = y2_init;

	for (int i = 0; i < N; i++) {

		t[i + 1] = t[i] + h;
		Yin_k1[0] = y1[i];
		Yin_k1[1] = y2[i];
		func(t[i], Yin_k1, K1); 
		double K1_y = K1[0]; // K1_y = dydt
		double K1_z = K1[1]; // K1_z = d2ydt2
		//printf("%f\n", K1_y); // dydt

		Yin_k2[0] = y1[i] + K1_y * h;
		Yin_k2[1] = y2[i] + K1_z * h;
		K2[0] = K1[0] + K1_y * h;
		K2[1] = K1[1] + K1_y * h;
		func(t[i] + h, Yin_k2, K2);
		double K2_y = K2[0];
		double K2_z = K2[1];

		y1[i + 1] = y1[i] + 0.5 * (K1_y + K2_y) * h;
		y2[i + 1] = y2[i] + 0.5 * (K1_z + K2_z) * h;
		//printf("%f", y1);
		
	}
}


```

**example**

A mass-spring damper (m-c-k) system is a second order ODEWhere  F(t) is the input force, y(t) is the displacement.

![img](file:///C:\Users\eunji\AppData\Local\Temp\msohtmlclip1\01\clip_image002.png)

![image-20221209004347705](C:\Users\eunji\AppData\Roaming\Typora\typora-user-images\image-20221209004347705.png)





## Linear Regression ## 

curve fitting : minimize the total residual 

```c++
int M = 6; // 6 개의 데이터 
double T_array[] = { 30, 40, 50, 60, 70, 80 };
double P_array[] = { 1.05, 1.07, 1.09, 1.14, 1.17, 1.21 };

Matrix T = arr2Mat(T_array, M, 1);
Matrix P = arr2Mat(P_array, M, 1);

Matrix z = linearRegression(T, P); // z = [a1, a0]

Matrix	linearRegression(Matrix x, Matrix y) {
	int mx = x.rows;
	int my = y.rows;
	int n = mx;

	double a1 = 0;
	double a0 = 0;

	double sx = 0; double sxx = 0; double sy = 0; double sxy = 0;
	double z_array[] = { a1,a0 };

	if (mx != my) {
		printf("\n****************************************************");
		printf("\n  ERROR!!: the length of matrix X and Y should be same");
		printf("\n****************************************************\n");
		return arr2Mat(z_array, 2, 1);
	}
	else {
		for (int k = 0; k < n; k++) {
			sx += x.at[k][0];
			sxx += x.at[k][0] * x.at[k][0];
			sxy += x.at[k][0] * y.at[k][0];
			sy += y.at[k][0];
		}

		double den = (n * sxx - sx * sx);
		a1 = (n * sxy - sx * sy) / den;
		a0 = (-sx * sxy + sxx * sy) / den;
	}
	z_array[0] = a1;
	z_array[1] = a0;
	return arr2Mat(z_array, 2, 1);
}
```



## Non-linear system with Jacobian 



```c++
void main() {
	
	
	double x0 = 3;
	double tol = 0.000001;
	double k = 0.0;
	double Nmax = 10;
	double ep = 1000;

	Matrix X1 = createMat(2, 1);
	Matrix H1 = createMat(2, 1);
	Matrix J1 = createMat(X1.rows, X1.rows);
	Matrix Jinv1 = createMat(X1.rows, X1.rows);
	Matrix F1 = createMat(X1.rows, X1.cols);

	X1.at[0][0] = 2.5;
	X1.at[1][0] = 2;
    
// Newton-Raphson Method Results


	do {
		/********************* Exercise 1 *********************/

		F1 = myFunc1(X1);
		J1 = myJacob1(X1);
		ep = vnorm(F1);

		//invMat(J, Jinv);
		//Jinv = mulConstMat(-1.0, Jinv);
		//H = multiplyMat(Jinv , F);
		Matrix L1 = createMat(J1.rows, J1.cols);
		Matrix U1 = createMat(J1.rows, J1.cols);
		LUdecomp(J1, L1, U1);
		Matrix F1_minus = mulConstMat(-1.0, F1);
		solveLU(L1, U1, F1_minus, H1);


		X1 = addMat(X1,H1);
		k++;

		/********************* Exercise 2 *********************/
		
		F2 = myFunc2(X2);
		J2 = myJacob2(X2);
		ep = vnorm(F2);

		//invMat(J, Jinv);
		//Jinv = mulConstMat(-1.0, Jinv);
		//H = multiplyMat(Jinv , F);
		Matrix L2 = createMat(J2.rows, J2.cols);
		Matrix U2 = createMat(J2.rows, J2.cols);
		LUdecomp(J2, L2, U2);
		Matrix F2_minus = mulConstMat(-1.0, F2);
		solveLU(L2, U2, F2_minus, H2);


		X2 = addMat(X2, H2);
		k++;



	} while (k < Nmax && ep > tol);
	printMat(X1, "X1\n"); printMat(X2, "X2\n");
}

F1 = myFunc1(X1);
		J1 = myJacob1(X1);
		ep = vnorm(F1);

		//invMat(J, Jinv);
		//Jinv = mulConstMat(-1.0, Jinv);
		//H = multiplyMat(Jinv , F);
		Matrix L1 = createMat(J1.rows, J1.cols);
		Matrix U1 = createMat(J1.rows, J1.cols);
		LUdecomp(J1, L1, U1);
		Matrix F1_minus = mulConstMat(-1.0, F1);
		solveLU(L1, U1, F1_minus, H1);


		X1 = addMat(X1,H1);
		k++;


Matrix myFunc1(Matrix X) {
	// matrix X is 2 by 1

	Matrix F = createMat(X.rows, X.cols);
	F.at[0][0] = X.at[1][0] - 0.5 * (exp(X.at[0][0] / 2) + exp(-X.at[0][0] / 2));
	F.at[1][0] = 9 * pow(X.at[0][0], 2) + 25 * pow(X.at[1][0], 2) - 225;
	return F;
}

Matrix myJacob1(Matrix X) {
	// matrix X is 2 by 2

	Matrix J = createMat(X.rows, X.rows);
	J.at[0][0] = -0.25 * exp(X.at[0][0] / 2 - exp(-X.at[0][0] / 2));
	J.at[0][1] = 1;
	J.at[1][0] = 18 * X.at[0][0];
	J.at[1][1] = 50 * X.at[1][0];
	return J;

}
/*==========================================================================*/
Matrix myFunc2(Matrix X) {
	// matrix X is 3 by 1, X = Z
	double x0 = 0;
	double y0 = 100;
	double x1 = 0;
	double y1 = -100;
	double th = X.at[0][0];
	double dx = X.at[1][0];
	double dy = X.at[2][0];

	double  x_new0 = 50;
	double y_new0 = 186.6025;
	double x_new1 = 150;
	double y_new1 = 13.3975;


	Matrix F = createMat(X.rows, X.cols);
	F.at[0][0] = x0 * cos(th) - y0 * sin(th) + dx - x_new0;
	F.at[1][0] = x0 * sin(th) + y0 * cos(th) + dy - y_new0;
	F.at[2][0] = x1 * cos(th) - y1 * sin(th) + dx - x_new1;
	return F;
}

Matrix myJacob2(Matrix X) {
	// matrix X is 3 by 3
	double x0 = 0;
	double y0 = 100;
	double x1 = 0;
	double y1 = -100;
	double th = X.at[0][0];

	Matrix J = createMat(X.rows, X.rows);
	J.at[0][0] = -x0 * sin(th) - y0 * cos(th);
	J.at[0][1] = 1;
	J.at[0][2] = 0;
	J.at[1][0] = x0 * cos(th) - y0 * sin(th);
	J.at[1][1] = 0;
	J.at[1][2] = 1;
	J.at[2][0] = -x1 * sin(th) - y1 * cos(th);
	J.at[2][1] = 1;
	J.at[2][2] = 0;

	return J;
}

```



## polyfit 

: curve fitting with a higher order poly

```c++
Matrix x = zeros(16, 1);
for (int i = 0; i < 16; i++) { x.at[i][0] = 0.4 * i;}

Matrix y = txt2Mat(path, "mat_y");
int n = 3;
//int n = 4;


Matrix z_poly = zeros(n + 1, 1);
z_poly = polyfit(x, y, z, n);

printMat(z, "z\n");
printMat(z_poly, "n= 3, z");
//printMat(z_poly, "n= 4, z");

Matrix	polyfit(Matrix x, Matrix y, int n) {

	Matrix z = zeros(n, 1);

	// n is order of polynominal 
	int mx = x.rows; // 16
	//int my = y.rows;
	//int m = sizeof(x) / sizeof(double);	printf("%d", m);

	Matrix sX = zeros(2 * n + 1, 1);
	Matrix sxy = zeros(n + 1, 1);
	Matrix s = zeros(n + 1, n + 1);
	Matrix sinv = zeros(n + 1, n + 1);
	Matrix b = zeros(n + 1, 1);
	double sxtemp = 0.0;
	double sxytemp = 0.0;

	for (int i = 0; i < 2 * n + 1; i++) {
		sxtemp = 0.0;
		for (int k = 0; k < mx; k++) {
			sxtemp = sxtemp + pow(x.at[k][0], i);
		}
		sX.at[i][0] = sxtemp;
	}
	for (int i = 0; i < n + 1; i++) {
		for (int j = 0; j < n + 1; j++) {
			s.at[i][j] = sX.at[2 * n - i - j][0];
		}
	}

	for (int j = n; j >= 0; j--) {
		sxytemp = 0.0;
		for (int k = 0; k < mx; k++) {
			sxytemp = sxytemp + y.at[k][0] * pow(x.at[k][0], j);
		}
		sxy.at[j][0] = sxytemp;
	}
	for (int i = 0; i < n + 1; i++) {
		b.at[i][0] = sxy.at[n - i][0];
	}

	invMat(s, sinv);
	z = multiplyMat(sinv, b);
	//printMat(z, "z");
	return z;
}

```











