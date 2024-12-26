======
Theory
======

QR decomposition
----------------

Square matrix
~~~~~~~~~~~~~

Let us consider a square matrix :math:`\mathbf{A}` of size :math:`n`, the QR decomposition of this matrix is:

.. math::
   \mathbf{A} = \mathbf{Q}\mathbf{R}
   :label: QR_decomposition_square

where:

 - :math:`\mathbf{Q}` is an orthogonal square matrix of size :math:`n` such as :math:`\mathbf{Q}^T\mathbf{Q} = \mathbf{Q}\mathbf{Q}^T = \mathbf{I}_n` with :math:`\mathbf{I}_n` the identity matrix of size :math:`n`;
 - :math:`\mathbf{R}` denotes an upper triangular square matrix of size :math:`n`.

The methods ``GetQRDecomposition`` and ``GetFullQRDecomposition`` return the matrices :math:`\mathbf{Q}` and :math:`\mathbf{R}`.

Rectangular matrix
~~~~~~~~~~~~~~~~~~

Let us consider a rectangular matrix :math:`A` of size :math:`m\times n`, the **full** decomposition of this matrix is:

.. math::
   \mathbf{A} = \mathbf{Q}\mathbf{R} = \mathbf{Q}\begin{pmatrix}\mathbf{R}_1 \\ \mathbf{0} \end{pmatrix}
   :label: QR_decomposition_rectangular_full

where:

 - :math:`\mathbf{Q}` is an orthogonal **square** matrix of size :math:`m` such as :math:`\mathbf{Q}^T\mathbf{Q} = \mathbf{Q}\mathbf{Q}^T = \mathbf{I}_n` with :math:`\mathbf{I}_m` the identity matrix of size :math:`n`;
 - :math:`\mathbf{R}` denotes an upper triangular **rectangular** matrix of size :math:`m\times n`;
 - :math:`\mathbf{R}_1` denotes an upper **square** triangular matrix of size :math:`n`;
 - :math:`\mathbf{0}` represents the null **rectangular** matrix of size :math:`(m-n)\times n`.

The matrix :math:`\mathbf{Q}` may be decomposed in two rectangular block matrices :math:`\mathbf{Q}_1` of size :math:`m\times n` and :math:`\mathbf{Q}_2` of size :math:`m\times (m - n)` such as:

.. math::
  \mathbf{Q} = \left(\mathbf{Q}_1, \mathbf{Q}_2\right)

The **thin** decomposition of :math:`\mathbf{A}` is:

.. math::
   \mathbf{A} = \mathbf{Q}_1\mathbf{R}_1
   :label: QR_decomposition_rectangular_thin

For a square matrix, these two decompositions are identical.

In **MathUtils**, the method ``GetQRDecomposition`` applies the thin decomposition and returns the matrices :math:`\mathbf{Q}_1` and :math:`\mathbf{R}_1` whereas the method ``GetFullQRDecomposition`` applies a full decomposition and returns the matrices :math:`\mathbf{Q}` and :math:`\mathbf{R}`.

.. _LS_problem:

Least-square problem
--------------------

Let us consider an overdetermined linear system :math:`\mathbf{A}\mathbf{x} = \mathbf{b}` with :math:`\mathbf{A}` a rectangular matrix of size :math:`m\times n` and :math:`\mathbf{b}` a vector of size :math:`m`. We assume the following condition:

.. math::
   m > n

The solving of this system is achieved in **MathUtils** by using a least-square method (``LeastSquareSolver``), based on the **Eigen** method: ``bdcSvd``:

.. math::
   \min\left\Vert\mathbf{A}\mathbf{x} - \mathbf{b}\right\Vert^2

The solution :math:`\mathbf{x}` satisfies:

.. math::
   \mathbf{x} = \mathbf{A}^{-1}\mathbf{b}
   :label: Solution_LS_problem

where :math:`\mathbf{A}^{-1}` is the pseudoinverse matrix of :math:`\mathbf{A}`.

Least-square problem subject to equality constraints
----------------------------------------------------

Let us consider the same problem as presented in :ref:`LS_problem` but we now add an equality constraint, so the problem becomes:

.. math::
   \min\left\Vert\mathbf{A}\mathbf{x} - \mathbf{b}\right\Vert^2 \: \text{subject to} \: \mathbf{C}\mathbf{x} = \mathbf{d}

with :math:`\mathbf{C}` a rectangular matrix of size :math:`p\times n` and :math:`\mathbf{d}` a vector of size :math:`p`.

For solving this problem, we follow the demonstration of [Grivet2016]_ in the section D.3.2. We start by doing a QR decomposition to the matrix :math:`\mathbf{C}^T` using :eq:`QR_decomposition_rectangular_full`:

.. math::
   \mathbf{C}^T = \mathbf{Q}\begin{pmatrix}\mathbf{R} \\ \mathbf{0} \end{pmatrix}

with :math:`\mathbf{Q}` a square matrix of size :math:`n` and :math:`\mathbf{R}` a square matrix of size :math:`p`.

Let us define the vector :math:`\mathbf{y}` of size :math:`p` and the vector :math:`\mathbf{z}` of size :math:`n-p` such as:

.. math::
   \mathbf{x} = \mathbf{Q}\begin{pmatrix} \mathbf{y} \\ \mathbf{z} \end{pmatrix}
   :label: Definition_y_z

By using the orthogonality of the matrix :math:`\mathbf{Q}`, we get:

.. math::
   \mathbf{C}\mathbf{x} = \mathbf{R}^T\mathbf{y}

Thus, the unknown :math:`\mathbf{y}` is determined from:

.. math::
   \mathbf{R}^T\mathbf{y} = \mathbf{d}
   :label: Linear_system_y

Let us define the two rectangular block matrices :math:`\mathbf{A}_1` of size :math:`m\times p` and :math:`\mathbf{A}_2` of size :math:`m\times (n - p)` such as:

.. math::
  \mathbf{A}\mathbf{Q} = \left(\mathbf{A}_1, \mathbf{A}_2\right)

We have:

.. math::
   \mathbf{A}\mathbf{x} = \mathbf{A}_1\mathbf{y} + \mathbf{A}_2\mathbf{z}

The unknown :math:`\mathbf{z}` is figured out by the following least-square problem where :math:`\mathbf{y}` is oibtained by :eq:`Linear_system_y`:

.. math::
   \min\left\Vert\mathbf{A}_2\mathbf{z} - \left(\mathbf{b} - \mathbf{A}_1\mathbf{y}\right)\right\Vert^2

Finally, the solution :math:`\mathbf{x}` is obtained from :eq:`Definition_y_z`.

This solution satisfies exactly the equality constraint :math:`\mathbf{C}\mathbf{x} = \mathbf{d}` but not exactly the linear system :math:`\mathbf{A}\mathbf{x} = \mathbf{b}`. Thus, contrairy to :ref:`LS_problem`, the solution does not satisfy :eq:`Solution_LS_problem`. Some errors appear because of the equality constraint [Grivet2016]_.

This method is used in ``LeastSquareSolverConstraint``.

Linear interpolation
--------------------

Let us consider a function :math:`f` over the segment :math:`[x_1;x_2]`. The values of :math:`f` at the ends of the segment are known. The linear interpolation of :math:`f` over the segment is:

.. math::
   f(x) = a_1 + a_2x

with:

.. math::
   \begin{cases}
       a_1 = \dfrac{x_2f(x_1) - x_1f(x_2)}{x_2 - x_1}\\
       a_2 = \dfrac{f(x_2) - f(x_1)}{x_2 - x_1}
   \end{cases}

This approach is used in the class ``Interp1d``.

Bilinear interpolation
----------------------

Let us consider a function :math:`f` over the set :math:`I = [x_1;x_2] \times [y_1;y_2]`. The values of :math:`f` at the ends of each segment are known. The bilinear interpolation of :math:`f` over :math:`I` is:

.. math::
   f(x, y) = a_1 + a_2x + a_3y + a_4xy

with:

.. math::
   \begin{cases}
       a_1 = \dfrac{x_2y_2f(x_1, y_1) - x_2y_1f(x_1, y_1) - x_1y_2f(x_2, y_1) + x_1y_1f(x_2, y_2)}{\Delta x\Delta y}\\
       a_2 = \dfrac{-y_2f(x_1, y_1) + y_1f(x_1, y_2) + y_2f(x_2, y_1) - y_1f(x_2, y_2)}{\Delta x\Delta y}\\
       a_3 = \dfrac{-x_2f(x_1, y_1) + x_1f(x_1, y_2) + x_2f(x_2, y_1) - x_1f(x_2, y_2)}{\Delta x\Delta y}\\
       a_4 = \dfrac{f(x_1, y_1) - f(x_1, y_2) - f(x_2, y_1) + f(x_2, y_2)}{\Delta x\Delta y}
   \end{cases}

where:

.. math::
   \begin{cases}
       \Delta x = x_2 - x_1\\
       \Delta y = y_2 - y_1
   \end{cases}

This approach is used in the class ``Interp2d``.

Legendre polynomials
--------------------

One definition of the Legendre polynomial of order :math:`n`, written :math:`P_n`, for :math:`x \in [-1, 1]` is:

.. math::
   \dfrac{d}{dx}\left[(1-x^2)\dfrac{d}{dx}P_n(x)\right] + n(n+1)P_n(x) = 0

This polynomial is real.

The Legendre polynomials obey the following recurrence relation:

.. math::
   \begin{cases}
      P_n(x) = \dfrac{1}{n}\left[(2n-1)xP_{n-1}(x) - (n-1)P_{n-2}(x)\right] \text{ for } n \geqslant 2\\
      P_0(x) = 1\\
      P_1(x) = x
   \end{cases}

The computation of the Legendre polynomials is achieved in the function ``Legendre_polynomial``.

Their differentiation is given by the following recurrence relation:

.. math::
   \begin{cases}
      P_n^{'}(x) = \dfrac{1}{x^2-1}\left[nxP_n(x) - nP_{n-1}(x)\right] \text{ for } n \geqslant 1\\
      P_0^{'}(x) = 0
   \end{cases}

The computation of the differentiation of the Legendre polynomials is achieved in the function ``Legendre_polynomial_derivative``.

Gamma function
--------------

The Euler's gamma function is defined by:

.. math::
   \Gamma(x) = \displaystyle \int_{0}^{+\infty} t^{x-1}e^{-t} dt

This function may be called with the function ``Gamma``.

For a positive interger :math:`n`, it comes:

.. math::
   \Gamma(n) = (n-1)!

The factorial may also be evaluated with the function ``Factorial``.

A particular value is:

.. math::
   \Gamma\left(n + \dfrac{1}{1}\right) = \dfrac{(2n)!}{2^{2n}n!}\sqrt{\pi}

For example:

.. math::
   \Gamma\left(\dfrac{1}{2}\right) &= \sqrt{\pi} \\ \Gamma\left(\dfrac{3}{2}\right) &= \dfrac{\sqrt{\pi}}{2} \\ \Gamma\left(\dfrac{5}{2}\right) &= \dfrac{3}{4}\sqrt{\pi}

Exponential integral :math:`Ei`
-------------------------------

The exponential integral :math:`Ei` is defined by [Abramowitz1964]_:

.. math::
   Ei(x) = \displaystyle -\int_{-x}^{+\infty}\dfrac{e^{-t}}{t} dt \text{ for } x > 0

This integral is evaluated in the function ``Ei``.

An approximation of order :math:`n` of this function was given by [Cody1969]_:

.. math::
   Ei(x) \approx \begin{cases} \ln\left(\dfrac{x}{x_0}\right) + (x - x_0)\dfrac{\displaystyle\sum_{j = 0}^n{}^{'} p_jT_j^{*}\left(\dfrac{x}{6}\right)}{\displaystyle\sum_{j = 0}^n{}^{'} q_jT_j^{*}\left(\dfrac{x}{6}\right)} \text{ for } 0 < x \leqslant 6 \\ \dfrac{e^x}{x}\left(\alpha_0 + \dfrac{\beta_0}{\alpha_1 + x + \dfrac{\beta_1}{\alpha_2 + x + \dfrac{\beta_2}{\alpha_3 + x + \ddots}}}\right)  \text{ for } \begin{cases} 6 < x \leqslant 12 \\ 12 < x \leqslant 24 \end{cases} \\ \dfrac{e^x}{x}\left[1 + \dfrac{1}{x}\left(\alpha_0 + \dfrac{\beta_0}{\alpha_1 + x + \dfrac{\beta_1}{\alpha_2 + x + \dfrac{\beta_2}{\alpha_3 + x + \ddots}}}\right)\right] \text{ for } x > 24 \end{cases}

where :math:`x_0` is the zero of :math:`Ei` and :math:`T_j^{*}` is a shifted Chebyshev polynomial  at the order :math:`j` (cf. :ref:`Chebyshev_polynomial`) defined by:

.. math::
   T_j^{*}(x) = T_j(2x-1)

.. note::
   The prime summation :math:`\sum{}^{'}` indicates only half of the first term is included.

The parameters :math:`x_0`, :math:`(p_j)_{0 \leqslant j \leqslant n}`, :math:`(q_j)_{0 \leqslant j \leqslant n}`, :math:`(\alpha_j)_{0 \leqslant j \leqslant n}` and :math:`(\beta_j)_{0 \leqslant j \leqslant n-1}` are given in [Cody1969]_. 

In **MathUtils**, the order is fixed to :math:`n = 9`.

The computation of this approximation of :math:`Ei` is achieved in the function ``Ei_approximation``.

The equivalent for :math:`x \to 0` is:

.. math::
   Ei(x) \sim \ln(x) + \gamma

The equivalent for :math:`x \to +\infty` is:

.. math::
   Ei(x) \sim \dfrac{e^{x}}{x}

It may be necessary to evaluate the quantity :math:`e^{-x}Ei(x)` for large :math:`x`, for example with the finite-depth Green's function for large water depth. This becomes impossible numerically as :math:`Ei` tends to infinity for large :math:`x`. Nevertheless, from the previous approximation, it comes:

.. math::
   e^{-x}Ei(x) \approx \dfrac{1}{x}\left[1 + \dfrac{1}{x}\left(\alpha_0 + \dfrac{\beta_0}{\alpha_1 + x + \dfrac{\beta_1}{\alpha_2 + x + \dfrac{\beta_2}{\alpha_3 + x + \ddots}}}\right)\right] \text{ for } x > 24

Which may be evaluated numerically without difficulty.

The computation of :math:`e^{-x}Ei(x)` is achieved in the function ``expEi``.

Exponential integral :math:`E_1`
--------------------------------

The exponential integral :math:`E_1` is defined by [Abramowitz1964]_:

.. math::
   E_1(x) = \displaystyle \int_{x}^{+\infty}\dfrac{e^{-t}}{t} dt \text{ for } x > 0

This integral is evaluated in the function ``E1``.

For :math:`x > 0`:

.. math::
   E_1(x) = -Ei(-x)

An approximation of order :math:`n` of this function was given by [Cody1968]_:

.. math::
   E_1(x) \approx \begin{cases} -\ln(x) + \dfrac{\displaystyle\sum_{j = 0}^n p_jx^j}{\displaystyle\sum_{j = 0}^n q_jx^j} \text{ for } 0 < x \leqslant 1 \\ e^{-x}\dfrac{\displaystyle\sum_{j = 0}^n p_jx^{-j}}{\displaystyle\sum_{j = 0}^n q_jx^{-j}}  \text{ for } 1 < x \leqslant 4 \\ \dfrac{e^{-x}}{x}\left[1 + \dfrac{1}{x}\dfrac{\displaystyle\sum_{j = 0}^n p_jx^{-j}}{\displaystyle\sum_{j = 0}^n q_jx^{-j}}\right] \text{ for } x > 4 \end{cases}

The coefficients :math:`(p_j)_{0 \leqslant j \leqslant n}` and :math:`(q_j)_{0 \leqslant j \leqslant n}` for each domain are given in [Cody1968]_. 

In **MathUtils**, the order is fixed to :math:`n = 6`, :math:`n = 8` and :math:`n = 9` for the three domains.

The computation of this approximation of :math:`E_1` is achieved in the function ``E1_approximation``.

The equivalent for :math:`x \to 0` is:

.. math::
   E_1(x) \sim -\ln(x) - \gamma

The equivalent for :math:`x \to +\infty` is:

.. math::
   E_1(x) \sim \dfrac{e^{-x}}{x}

Struve functions
----------------

The Struve function of order :math:`n`, :math:`H_n`, is defined in [Abramowitz1964]_ (chapter 12). Efficient approximations for the orders zero and one are provided by [Newman1984]_.

The zero-order Struve function is evaluted using:

.. math::
   H_0(x) \approx \begin{cases}
      \displaystyle \sum_{j = 1}^6 a_j\left(\dfrac{x}{3}\right)^{2j-1} \text{ for } 0 \leqslant x \leqslant 3\\
      Y_0(x) + \dfrac{2}{\pi x}\dfrac{\displaystyle \sum_{j = 0}^3 b_j\left(\dfrac{3}{x}\right)^{2j}}{\displaystyle\sum_{j = 0}^3 c_j\left(\dfrac{3}{x}\right)^{2j}} \text{ for } x > 3
   \end{cases}

with:

.. math::
   \begin{cases}
      a_1 = 1.909859164\\
      a_2 = -1.909855001\\
      a_3 = 0.687514637\\
      a_4 = -0.126164557\\
      a_5 = 0.013828813\\
      a_6 = -0.000876918
   \end{cases}

and

.. math::
   \begin{cases}
      b_0 = 0.99999906\\
      b_1 = 4.77228920\\
      b_2 = 3.85542044\\
      b_3 = 0.32303607\\
      c_0 = 1\\
      c_1 = 4.88331068\\
      c_2 = 4.28957333\\
      c_3 = 0.52120508
   \end{cases}

For :math:`x \to 0`, there is:

.. math::
   H_0(x) \sim \dfrac{2}{\pi}x

The first-order Struve function is evaluted using:

.. math::
   H_1(x) \approx \begin{cases}
      \displaystyle \sum_{j = 1}^6 d_j\left(\dfrac{x}{3}\right)^{2j} \text{ for } 0 \leqslant x \leqslant 3\\
      Y_1(x) + \dfrac{2}{\pi}\dfrac{\displaystyle \sum_{j = 0}^3 e_j\left(\dfrac{3}{x}\right)^{2j}}{\displaystyle\sum_{j = 0}^3 f_j\left(\dfrac{3}{x}\right)^{2j}} \text{ for } x > 3
   \end{cases}

with:

.. math::
   \begin{cases}
      d_1 = 1.909859286\\
      d_2 = -1.145914713\\
      d_3 = 0.294656958\\
      d_4 = -0.042070508\\
      d_5 = 0.003785727\\
      d_6 = -0.000207183
   \end{cases}

and

.. math::
   \begin{cases}
      e_0 = 1.00000004\\
      e_1 = 3.92205313\\
      e_2 = 2.64893033\\
      e_3 = 0.27450895\\
      f_0 = 1\\
      f_1 = 3.81095112\\
      f_2 = 2.26216956\\
      f_3 = 0.10885141
   \end{cases}

:math:`Y_0` and :math:`Y_1` represent the Bessel functions of second kind of order zero and one.

For :math:`x \to 0`, there is:

.. math::
   H_1(x) \sim \dfrac{2}{3\pi}x^2

The derivative of :math:`H_0` is given by [Abramowitz1964]_:

.. math::
   H_0^{'}(x) = \dfrac{2}{\pi} - H_1(x)

The computation of the zero-order and first-order of the Struve functions are achieved in the functions ``Struve_zero_order`` and ``Struve_first_order``. Regarding the derivative of :math:`H_0`, the function to use is ``Struve_zero_order_derivative``.

For :math:`x \to 0`, there is:

.. math::
   H_0^{'}(x) \sim \dfrac{2}{\pi}

.. _Chebyshev_polynomial:

Chebyshev polynomials
---------------------

The Chebyshev polynomials, written :math:`T_n` at the order :math:`n`, are defind by the following recurrence relation for :math:`x \in [-1, 1]`:

.. math::
   \begin{cases}
      T_{n+2}(x) = 2xT_{n+1}(x) - T_{n}(x) \text{ for } n \geqslant 0\\
      T_0(x) = 1\\
      T_1(x) = x
   \end{cases}

The zeros of :math:`T_{n+1}` for :math:`x \in [-1, 1]` and :math:`i \in [0, n]` are:

.. math::
   x_i = \cos\left[\dfrac{\pi}{2}\left(\dfrac{2i+1}{n+1}\right)\right]


If :math:`x \in [a, b]`, then, by affine transformation:

.. math::
   x \rightarrow \dfrac{2}{b-a}\left(x - \dfrac{b+a}{2}\right)

The zeros of :math:`T_{n+1}` for :math:`x \in [a, b]` and :math:`i \in [0, n]` are:

.. math::
   x_i = \left(\dfrac{b-a}{2}\right)\cos\left[\dfrac{\pi}{2}\left(\dfrac{2i+1}{n+1}\right)\right] + \dfrac{b+a}{2}

The computation of the Chebyshev polynomials is achieved in the functions ``Chebyshev_polynomial`` and ``Chebyshev_polynomial_next``.

Single Chebyshev series approximation
-------------------------------------

The single Chebyshev series approximation of order :math:`m` of the function :math:`f` defined over :math:`[x_{min}, x_{max}]` is expressed by [Basu1973]_:

.. math::
   \displaystyle f(x) \approx \sum_{i = 0}^m a_{i}T_{i}(\tilde{x})
   :label: Single_Chebyshev_approx

with:

.. math::
   \begin{cases} \tilde{x} = \dfrac{2}{x_{max} - x_{min}}\left[x - \left(\dfrac{x_{max} + x_{min}}{2}\right)\right] \in [-1, 1] \\ a_{i} = \begin{cases} \displaystyle \dfrac{1}{m+1}\sum_{r = 0}^m f(x_r)T_{i}(\tilde{x}_r) \text{ if } i = 0 \\ \displaystyle \dfrac{2}{m+1}\sum_{r = 0}^m f(x_r)T_{i}(\tilde{x}_r) \text{ if } i \neq 0 \end{cases} \\ \tilde{x}_r = \cos\left[\dfrac{\pi}{2}\left(\dfrac{2r+1}{m+1}\right)\right] \in [-1, 1] \\ x_r = \left(\dfrac{x_{max}-x_{min}}{2}\right)\tilde{x}_r + \dfrac{x_{max}+x_{min}}{2} \in [x_{min}, x_{max}]
   \end{cases}

:math:`T_i` represents the Chebyshev polynomial of order :math:`i`.

The coefficients :math:`(a_{i})_{0 \leqslant i \leqslant m}` must be computed in a first step before evaluating :eq:`Single_Chebyshev_approx` for any value of :math:`x`. 

At the points :math:`(x_r)_{0 \leqslant r \leqslant m}`, by definition, it yields:

.. math::
   \displaystyle f(x_r) = \sum_{i = 0}^m a_{i}T_{i}(\tilde{x}_r)

Regarding the derivative, it comes:

.. math::
   \displaystyle \dfrac{\partial f}{\partial x}(x) \approx \left(\dfrac{2}{x_{max} - x_{min}}\right)\sum_{i = 0}^m a_{i}\dfrac{\partial T_i}{\partial x}(\tilde{x})

.. note::
   If the function :math:`f` is defined over a half-open line segment :math:`[x_{min}, +\infty[`, then the following modifications are necessary [Chen1993]_:

      .. math::
          \begin{cases}
              \tilde{x} = 1 - 2\dfrac{x_{min}}{x} \in [-1, 1] \\ x_r = 2\dfrac{x_{min}}{1 - \tilde{x}_r} \in [x_{min}, +\infty[ \\ \displaystyle \dfrac{\partial f}{\partial x}(x) \approx \dfrac{2x_{min}}{x^2}\sum_{i = 0}^m a_{i}\dfrac{\partial T_i}{\partial x}(\tilde{x})
          \end{cases}

The single Chebyshev series approximation is performed with the base class ``ChebyshevSeries1dBase`` and its derived classes.

Double Chebyshev series approximation
-------------------------------------

The double Chebyshev series approximation of order :math:`m\times n` of the function :math:`f` defined over :math:`[x_{min}, x_{max}]\times[y_{min}, y_{max}]` is expressed by [Basu1973]_:

.. math::
   \displaystyle f(x,y) \approx \sum_{i = 0}^m\sum_{j = 0}^n a_{ij}T_{i,j}(\tilde{x}, \tilde{y})
   :label: Double_Chebyshev_approx

with:

.. math::
   \begin{cases} \tilde{x} = \dfrac{2}{x_{max} - x_{min}}\left[x - \left(\dfrac{x_{max} + x_{min}}{2}\right)\right] \in [-1, 1] \\ \tilde{y} = \dfrac{2}{y_{max} - y_{min}}\left[y - \left(\dfrac{y_{max} + y_{min}}{2}\right)\right] \in [-1, 1] \\ a_{ij} = \begin{cases} \displaystyle \dfrac{1}{(m+1)(n+1)}\sum_{r = 0}^m\sum_{s = 0}^nf(x_r, y_s)T_{i,j}(\tilde{x}_r, \tilde{y}_s) \text{ if } \begin{cases} i = 0 \\ j = 0 \end{cases} \\ \displaystyle \dfrac{2}{(m+1)(n+1)}\sum_{r = 0}^m\sum_{s = 0}^nf(x_r, y_s)T_{i,j}(\tilde{x}_r, \tilde{y}_s) \text{ if } \begin{cases} i = 0 \\ j \neq 0 \end{cases} \text{ or } \begin{cases} i \neq 0 \\ j = 0 \end{cases} \\ \displaystyle \dfrac{4}{(m+1)(n+1)}\sum_{r = 0}^m\sum_{s = 0}^nf(x_r, y_s)T_{i,j}(\tilde{x}_r, \tilde{y}_s) \text{ otherwise }\end{cases} \\ \tilde{x}_r = \cos\left[\dfrac{\pi}{2}\left(\dfrac{2r+1}{m+1}\right)\right] \in [-1, 1] \\ \tilde{y}_s = \cos\left[\dfrac{\pi}{2}\left(\dfrac{2s+1}{n+1}\right)\right] \in [-1, 1] \\ x_r = \left(\dfrac{x_{max}-x_{min}}{2}\right)\tilde{x}_r + \dfrac{x_{max}+x_{min}}{2} \in [x_{min}, x_{max}] \\ y_s = \left(\dfrac{y_{max}-y_{min}}{2}\right)\tilde{y}_s + \dfrac{y_{max}+y_{min}}{2} \in [y_{min}, y_{max}] \\ T_{i,j}(\tilde{x}, \tilde{y}) = T_i(\tilde{x})T_j(\tilde{y})\end{cases}

:math:`T_i` represents the Chebyshev polynomial of order :math:`i`.

The coefficients :math:`(a_{ij})_{0 \leqslant i \leqslant m \\ 0 \leqslant j \leqslant n}` must be computed in a first step before evaluating :eq:`Double_Chebyshev_approx` for any value of :math:`x` and :math:`y`. 

At the points :math:`(x_r, y_s)_{0 \leqslant r \leqslant m \\ 0 \leqslant s \leqslant n}`, by definition, it yields:

.. math::
   \displaystyle f(x_r, y_s) = \sum_{i = 0}^m\sum_{j = 0}^n a_{ij}T_{i,j}(\tilde{x}_r, \tilde{y}_s)

Regarding the partial derivatives, it comes:

.. math::
   \begin{cases}
   \displaystyle \dfrac{\partial f}{\partial x}(x,y) \approx \left(\dfrac{2}{x_{max} - x_{min}}\right)\sum_{i = 0}^m\sum_{j = 0}^n a_{ij}\dfrac{\partial T_i}{\partial x}(\tilde{x})T_j(\tilde{y})\\
   \displaystyle \dfrac{\partial f}{\partial y}(x,y) \approx \left(\dfrac{2}{y_{max} - y_{min}}\right)\sum_{i = 0}^m\sum_{j = 0}^n a_{ij}T_i(\tilde{x})\dfrac{\partial T_j}{\partial y}(\tilde{y})
   \end{cases}

.. note::
   If the function :math:`f` is defined over half-open line segments :math:`[x_{min}, +\infty[\times[y_{min}, +\infty[`, then the following modifications are necessary [Chen1993]_:

      .. math::
          \begin{cases}
              \tilde{x} = 1 - 2\dfrac{x_{min}}{x} \in [-1, 1] \\ \tilde{y} = 1 - 2\dfrac{y_{min}}{y} \in [-1, 1] \\ x_r = 2\dfrac{x_{min}}{1 - \tilde{x}_r} \in [x_{min}, +\infty[ \\ y_s = 2\dfrac{y_{min}}{1 - \tilde{y}_s} \in [y_{min}, +\infty[\\ \displaystyle \dfrac{\partial f}{\partial x}(x,y) \approx \dfrac{2x_{min}}{x^2}\sum_{i = 0}^m\sum_{j = 0}^n a_{ij}\dfrac{\partial T_i}{\partial x}(\tilde{x})T_j(\tilde{y})\\ \displaystyle \dfrac{\partial f}{\partial y}(x,y) \approx \dfrac{2y_{min}}{y^2}\sum_{i = 0}^m\sum_{j = 0}^n a_{ij}T_i(\tilde{x})\dfrac{\partial T_j}{\partial y}(\tilde{y})
          \end{cases}

The double Chebyshev series approximation is performed with the base class ``ChebyshevSeries2dBase`` and its derived classes.

Triple Chebyshev series approximation
-------------------------------------

The triple Chebyshev series approximation of order :math:`m\times n\times p` of the function :math:`f` defined over :math:`[x_{min}, x_{max}]\times[y_{min}, y_{max}]\times[z_{min}, z_{max}]` is obtained from [Mackay2019]_ and the generalization of the previous section:

.. math::
   \displaystyle f(x,y, z) \approx \sum_{i = 0}^m\sum_{j = 0}^n\sum_{k = 0}^p a_{ijk}T_{i,j,k}(\tilde{x}, \tilde{y}, \tilde{z})
   :label: Triple_Chebyshev_approx

with:

.. math::
   \begin{cases} \tilde{x} = \dfrac{2}{x_{max} - x_{min}}\left[x - \left(\dfrac{x_{max} + x_{min}}{2}\right)\right] \in [-1, 1] \\ \tilde{y} = \dfrac{2}{y_{max} - y_{min}}\left[y - \left(\dfrac{y_{max} + y_{min}}{2}\right)\right] \in [-1, 1] \\ \tilde{z} = \dfrac{2}{z_{max} - z_{min}}\left[z - \left(\dfrac{z_{max} + z_{min}}{2}\right)\right] \in [-1, 1] \\ a_{ijk} = \begin{cases} \displaystyle \dfrac{1}{(m+1)(n+1)(p+1)}\sum_{r = 0}^m\sum_{s = 0}^n\sum_{t = 0}^pf(x_r, y_s, z_t)T_{i,j,k}(\tilde{x}_r, \tilde{y}_s, \tilde{z}_t) \text{ if } \begin{cases} i = 0 \\ j = 0 \\ k = 0 \end{cases} \\ \displaystyle \dfrac{2}{(m+1)(n+1)(p+1)}\sum_{r = 0}^m\sum_{s = 0}^n\sum_{t = 0}^pf(x_r, y_s, z_t)T_{i,j,k}(\tilde{x}_r, \tilde{y}_s, \tilde{z}_t) \text{ if } \begin{cases} i = 0 \\ j \neq 0 \\ k = 0 \end{cases} \text{ or } \begin{cases} i \neq 0 \\ j = 0 \\ k = 0 \end{cases}  \text{ or } \begin{cases} i = 0 \\ j = 0 \\ k \neq 0 \end{cases} \\ \displaystyle \dfrac{4}{(m+1)(n+1)(p+1)}\sum_{r = 0}^m\sum_{s = 0}^n\sum_{t = 0}^pf(x_r, y_s, z_t)T_{i,j,k}(\tilde{x}_r, \tilde{y}_s, \tilde{z}_t) \text{ if } \begin{cases} i = 0 \\ j \neq 0 \\ k \neq 0 \end{cases} \text{ or } \begin{cases} i \neq 0 \\ j = 0 \\ k \neq 0 \end{cases} \text{ or } \begin{cases} i \neq 0 \\ j \neq 0 \\ k = 0 \end{cases} \\ \displaystyle \dfrac{8}{(m+1)(n+1)(p+1)}\sum_{r = 0}^m\sum_{s = 0}^n\sum_{t = 0}^pf(x_r, y_s, z_t)T_{i,j,k}(\tilde{x}_r, \tilde{y}_s, \tilde{z}_t) \text{ otherwise }\end{cases} \\ \tilde{x}_r = \cos\left[\dfrac{\pi}{2}\left(\dfrac{2r+1}{m+1}\right)\right] \in [-1, 1] \\ \tilde{y}_s = \cos\left[\dfrac{\pi}{2}\left(\dfrac{2s+1}{n+1}\right)\right] \in [-1, 1] \\ \tilde{z}_t = \cos\left[\dfrac{\pi}{2}\left(\dfrac{2t+1}{p+1}\right)\right] \in [-1, 1] \\ x_r = \left(\dfrac{x_{max}-x_{min}}{2}\right)\tilde{x}_r + \dfrac{x_{max}+x_{min}}{2} \in [x_{min}, x_{max}] \\ y_s = \left(\dfrac{y_{max}-y_{min}}{2}\right)\tilde{y}_s + \dfrac{y_{max}+y_{min}}{2} \in [y_{min}, y_{max}] \\ z_t = \left(\dfrac{z_{max}-z_{min}}{2}\right)\tilde{z}_t + \dfrac{z_{max}+z_{min}}{2} \in [z_{min}, z_{max}] \\ T_{i,j,k}(\tilde{x}, \tilde{y},\tilde{z}) = T_i(\tilde{x})T_j(\tilde{y})T_k(\tilde{z})\end{cases}

:math:`T_i` represents the Chebyshev polynomial of order :math:`i`.

The coefficients :math:`(a_{ijk})_{0 \leqslant i \leqslant m \\ 0 \leqslant j \leqslant n \\ 0 \leqslant k \leqslant p}` must be computed in a first step before evaluating :eq:`Triple_Chebyshev_approx` for any value of :math:`x`, :math:`y` and :math:`z`. 

At the points :math:`(x_r, y_s, z_t)_{0 \leqslant r \leqslant m \\ 0 \leqslant s \leqslant n \\ 0 \leqslant t \leqslant p}`, by definition, it yields:

.. math::
   \displaystyle f(x_r, y_s, z_t) = \sum_{i = 0}^m\sum_{j = 0}^n\sum_{k = 0}^p a_{ijk}T_{i,j,k}(\tilde{x}_r, \tilde{y}_s, \tilde{z}_t)

Regarding the partial derivatives, it comes:

.. math::
   \begin{cases}
   \displaystyle \dfrac{\partial f}{\partial x}(x,y,z) \approx \left(\dfrac{2}{x_{max} - x_{min}}\right)\sum_{i = 0}^m\sum_{j = 0}^n\sum_{k = 0}^p a_{ijk}\dfrac{\partial T_i}{\partial x}(\tilde{x})T_j(\tilde{y})T_k(\tilde{z})\\
   \displaystyle \dfrac{\partial f}{\partial y}(x,y,z) \approx \left(\dfrac{2}{y_{max} - y_{min}}\right)\sum_{i = 0}^m\sum_{j = 0}^n\sum_{k = 0}^p a_{ijk}T_i(\tilde{x})\dfrac{\partial T_j}{\partial y}(\tilde{y})T_k(\tilde{z})\\
   \displaystyle \dfrac{\partial f}{\partial z}(x,y,z) \approx \left(\dfrac{2}{z_{max} - z_{min}}\right)\sum_{i = 0}^m\sum_{j = 0}^n\sum_{k = 0}^p a_{ijk}T_i(\tilde{x})T_j(\tilde{y})\dfrac{\partial T_k}{\partial z}(\tilde{z})
   \end{cases}

.. note::
   If the function :math:`f` is defined over half-open line segments :math:`[x_{min}, +\infty[\times[y_{min}, +\infty[\times[z_{min}, +\infty[`, then the following modifications are necessary [Chen1993]_:

      .. math::
          \begin{cases}
              \tilde{x} = 1 - 2\dfrac{x_{min}}{x} \in [-1, 1] \\ \tilde{y} = 1 - 2\dfrac{y_{min}}{y} \in [-1, 1] \\ \tilde{z} = 1 - 2\dfrac{z_{min}}{z} \in [-1, 1] \\ x_r = 2\dfrac{x_{min}}{1 - \tilde{x}_r} \in [x_{min}, +\infty[ \\ y_s = 2\dfrac{y_{min}}{1 - \tilde{y}_s} \in [y_{min}, +\infty[ \\ z_t = 2\dfrac{z_{min}}{1 - \tilde{z}_t} \in [z_{min}, +\infty[ \\ \displaystyle \dfrac{\partial f}{\partial x}(x,y,z) \approx \dfrac{2x_{min}}{x^2}\sum_{i = 0}^m\sum_{j = 0}^n\sum_{k = 0}^p a_{ijk}\dfrac{\partial T_i}{\partial x}(\tilde{x})T_j(\tilde{y})T_k(\tilde{z})\\ \displaystyle \dfrac{\partial f}{\partial y}(x,y,z) \approx \dfrac{2y_{min}}{y^2}\sum_{i = 0}^m\sum_{j = 0}^n\sum_{k = 0}^p a_{ijk}T_i(\tilde{x})\dfrac{\partial T_j}{\partial y}(\tilde{y})T_k(\tilde{z})\\ \displaystyle \dfrac{\partial f}{\partial z}(x,y,z) \approx \dfrac{2z_{min}}{z^2}\sum_{i = 0}^m\sum_{j = 0}^n\sum_{k = 0}^p a_{ijk}T_i(\tilde{x})T_j(\tilde{y})\dfrac{\partial T_k}{\partial z}(\tilde{z})
          \end{cases}

The triple Chebyshev series approximation is performed with the base class ``ChebyshevSeries3dBase`` and its derived classes.

Horner's method
---------------

Let us define a polynomial :math:`P` of order :math:`n` such as:

.. math::
   \displaystyle P(x) = \sum_{k = 0}^n a_kx^k


The optimal algorithm for polynomial evaluation is the Horner's method. It requires the following rearrangement:

.. math::
   \displaystyle P(x) = a_0 + x(a_1 + x(\dots + x(a_{n-1} + xa_n)))

Starting from last coefficients, only :math:`n` multiplications and :math:`n` additions are required.

Regading the derivative of :math:`P`:

.. math::
   \displaystyle P'(x) = \sum_{k = 1}^n ka_kx^{k-1} = \sum_{k = 0}^{n-1} (k+1)a_{k+1}x^k

The method is applied using:

.. math::
   \begin{cases}
      \displaystyle P'(x) = \sum_{k = 0}^{n-1} b_kx^k\\
      b_k = (k+1)a_{k+1}
   \end{cases}

The functions ``Horner`` and ``Horner_derivative`` apply this method.

Conversion of Chebyshev series into power series
------------------------------------------------

It is interesting to convert Chebyshev series into power series for using the Horner's method and evaluating derivatives. Let us consider a function :math:`f` which is approximated by a double Chebyshev series of order :math:`m \times n` over :math:`[x_{min}, x_{max}] \times [y_{min}, y_{max}]` and converted into a power series:

.. math::
   \displaystyle f(x,y) \approx \sum_{i = 0}^m\sum_{j = 0}^n a_{ij}T_{i,j}(\tilde{x}, \tilde{y}) \approx \sum_{i = 0}^m\sum_{j = 0}^n b_{ij}\tilde{x}^i\tilde{y}^j

with:

.. math::
   \begin{cases} \tilde{x} = \dfrac{2}{x_{max} - x_{min}}\left[x - \left(\dfrac{x_{max} + x_{min}}{2}\right)\right] \in [-1, 1] \\ \tilde{y} = \dfrac{2}{y_{max} - y_{min}}\left[y - \left(\dfrac{y_{max} + y_{min}}{2}\right)\right] \in [-1, 1] \end{cases}

The series are evaluated using the Horner's method.

The coefficients :math:`(b_{ij})_{0 \leqslant i \leqslant m \\ 0 \leqslant j \leqslant n}` are expressed by [Chen1993]_:

.. math::
   \displaystyle b_{ij} = \sum_{r = i}^m\sum_{s = j}^n a_{rs}\lambda_{ri}\lambda_{sj}

with:

.. math::
   \lambda_{ri} = \begin{cases} 1 \text{ if } r = i = 0 \\ 0 \text{ if } (r + i) \text{ is odd} \\ \displaystyle (-1)^{\left[\dfrac{r - i}{2}\right]} 2^{i-1}\left(\dfrac{r\left[\dfrac{r + i}{2}-1\right]!}{\left[\dfrac{r - i}{2}\right]!i!}\right) \text{ otherwise}\end{cases}

where :math:`[x]` represents the floor function.

.. note::
   If the function :math:`f` is defined over half-open line segments :math:`[x_{min}, +\infty[\times[y_{min}, +\infty[`, then the following modifications are necessary:

      .. math::
          \begin{cases}
              \tilde{x} = 1 - 2\dfrac{x_{min}}{x} \in [-1, 1] \\ \tilde{y} = 1 - 2\dfrac{y_{min}}{y} \in [-1, 1]
          \end{cases}

.. warning::
   This algorithm presents numerical inaccuracies for large orders. Consequently, each order must be **lower or equal to 18** [Boyd2002]_.

If the function :math:`f` is approximated by a single Chebyshev series, then:

.. math::
   \displaystyle f(x) \approx \sum_{i = 0}^m a_{i}T_{i}(\tilde{x}) \approx \sum_{i = 0}^m b_{i}\tilde{x}^i

with:

.. math::
   \displaystyle b_{i} = \sum_{r = i}^m a_{r}\lambda_{ri}

If the function :math:`f` is approximated by a triple Chebyshev series, then:

.. math::
   \displaystyle f(x,y,z) \approx \sum_{i = 0}^m\sum_{j = 0}^n\sum_{k = 0}^p a_{ijk}T_{i,j,k}(\tilde{x}, \tilde{y}, \tilde{z}) \approx \sum_{i = 0}^m\sum_{j = 0}^n\sum_{k = 0}^p b_{ijk}\tilde{x}^i\tilde{y}^j\tilde{z}^k

with:

.. math::
   \displaystyle b_{ijk} = \sum_{r = i}^m\sum_{s = j}^n\sum_{t = k}^p a_{rst}\lambda_{ri}\lambda_{sj}\lambda_{tk}

These conversions are performed with the base classes ``PowerSeries1dBase``, ``PowerSeries2dBase`` and ``PowerSeries3dBase`` and their derived classes.

.. note::
   Let us consider a triple power series with a constant value of :math:`z`, the summation over :math:`k` may be made prior to the evaluations of :math:`f` for any values of :math:`x` and :math:`y`, such as:
   
   .. math::
      \displaystyle f(x,y,z) \approx \sum_{i = 0}^m\sum_{j = 0}^n\sum_{k = 0}^p b_{ijk}\tilde{x}^i\tilde{y}^j\tilde{z}^k \approx \sum_{i = 0}^m\sum_{j = 0}^n c_{ij}\tilde{x}^i\tilde{y}^j
      
   with:
   
   .. math::
      c_{ij} = \sum_{k = 0}^p b_{ijk}\tilde{z}^k
      
   In this case, triple power series become double power series.

.. [Abramowitz1964] M. Abramowitz and I. A. Stegun. Handbook of Mathematical functions with formulas, graphs and mathematical tables. Government Printing Office, Washington and Dover, New York, 1964.

.. [Cody1968] Cody W. J. and Thacher H. C. Rational Chebyshev Approximations for the Exponential Integral :math:`E_1(x)`. Mathematics of Computation, 22(103):641-649, 1968.

.. [Cody1969] Cody W. J. and Thacher H. C. Chebyshev Approximations for the Exponential Integral :math:`Ei(x)`. Mathematics of Computation, 23(106):289-303, 1969.

.. [Basu1973] N. K. Basu. On double Chebyshev series approximation. SIAM Journal on Numerical Analysis, 10(3):493-505, 1973.

.. [Newman1984] J. N. Newman. Approximations for the Bessel and Struve functions. Mathematics of Computation, 43(168):551-556, 1984.

.. [Chen1993] X. Chen. Evaluation de la fonction de Green du problème de diffraction / radiation en profondeur d’eau finie. Proceedings of the 4ème Journées de l’Hydrodynamique (JH1993), Nantes, France, 1993.

.. [Boyd2002] J. P. Boyd. Computing zeros on a real interval through Chebyshev expansion and polynomial rootfinding. SIAM Journal on Numerical Analysis, 40(5):1666-1682, 2002.

.. [Grivet2016] S. Grivet-Talocia and B. Gustavsen. Passive macromodeling. Theory and applications. 2016.

.. [Mackay2019] E. Mackay. Consistent expressions for the free-surface Green function in finite water depth. Applied Ocean Research, 93, 2019.
