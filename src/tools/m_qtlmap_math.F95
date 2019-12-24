MODULE LogBeta
! Log of the beta function
! Includes log of the gamma function

IMPLICIT NONE


CONTAINS


FUNCTION betaln(a0, b0) RESULT(fn_val)
!-----------------------------------------------------------------------
!     EVALUATION OF THE LOGARITHM OF THE BETA FUNCTION
!-----------------------------------------------------------------------
!     E = 0.5*LN(2*PI)
!--------------------------
IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(15, 100)

REAL (dp), INTENT(IN) :: a0, b0
REAL (dp)             :: fn_val

! Local variables
REAL (dp), PARAMETER  :: e = .918938533204673_dp
REAL (dp)             :: a, b, c, h, u, v, w, z
INTEGER               :: i, n
REAL (dp), PARAMETER  :: half = 0.5_dp, one = 1.0_dp, two = 2.0_dp,  &
                         eight = 8.0_dp
!--------------------------
a = MIN(a0, b0)
b = MAX(a0, b0)
IF (a < eight) THEN
  IF (a < one) THEN
!-----------------------------------------------------------------------
!                   PROCEDURE WHEN A < 1
!-----------------------------------------------------------------------
    IF (b < eight) THEN
      fn_val = gamln(a) + (gamln(b) - gamln(a+b))
      RETURN
    END IF
    fn_val = gamln(a) + algdiv(a,b)
    RETURN
  END IF
!-----------------------------------------------------------------------
!                PROCEDURE WHEN 1 <= A < 8
!-----------------------------------------------------------------------
  IF (a <= two) THEN
    IF (b <= two) THEN
      fn_val = gamln(a) + gamln(b) - gsumln(a,b)
      RETURN
    END IF
    w = 0.0
    IF (b < eight) GO TO 20
    fn_val = gamln(a) + algdiv(a,b)
    RETURN
  END IF

!                REDUCTION OF A WHEN B <= 1000

  IF (b > 1000.0_dp) GO TO 40
  n = a - one
  w = one
  DO i = 1, n
    a = a - one
    h = a / b
    w = w * (h/(one + h))
  END DO
  w = LOG(w)
  IF (b >= eight) THEN
    fn_val = w + gamln(a) + algdiv(a,b)
    RETURN
  END IF

!                 REDUCTION OF B WHEN B < 8

  20 n = b - one
  z = one
  DO i = 1, n
    b = b - one
    z = z * (b/(a+b))
  END DO
  fn_val = w + LOG(z) + (gamln(a) + (gamln(b) - gsumln(a,b)))
  RETURN

!                REDUCTION OF A WHEN B > 1000

  40 n = a - one
  w = one
  DO i = 1, n
    a = a - one
    w = w * (a/(one+a/b))
  END DO
  fn_val = (LOG(w) - n*LOG(b)) + (gamln(a) + algdiv(a,b))
  RETURN
END IF
!-----------------------------------------------------------------------
!                   PROCEDURE WHEN A >= 8
!-----------------------------------------------------------------------
w = bcorr(a,b)
h = a / b
c = h / (one + h)
u = -(a-half) * LOG(c)
v = b * alnrel(h)
IF (u > v) THEN
  fn_val = (((-half*LOG(b) + e) + w) - v) - u
  RETURN
END IF
fn_val = (((-half*LOG(b) + e) + w) - u) - v
RETURN
END FUNCTION betaln



FUNCTION gamln(a) RESULT(fn_val)
!-----------------------------------------------------------------------
!            EVALUATION OF LN(GAMMA(A)) FOR POSITIVE A
!-----------------------------------------------------------------------
!     WRITTEN BY ALFRED H. MORRIS
!          NAVAL SURFACE WARFARE CENTER
!          DAHLGREN, VIRGINIA
!--------------------------
!     D = 0.5*(LN(2*PI) - 1)
!--------------------------
IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(15, 100)

REAL (dp), INTENT(IN) :: a
REAL (dp)             :: fn_val

! Local variables
REAL (dp), PARAMETER :: d = 0.418938533204673_dp, c0 = 0.833333333333333D-01,   &
              c1 = -0.277777777760991D-02, c2 = 0.793650666825390D-03,   &
              c3 = -0.595202931351870D-03, c4 = 0.837308034031215D-03,   &
              c5 = -0.165322962780713D-02
REAL (dp)  :: t, w
INTEGER    :: i, n
!-----------------------------------------------------------------------
IF (a <= 0.8_dp) THEN
  fn_val = gamln1(a) - LOG(a)
  RETURN
END IF
IF (a <= 2.25_dp) THEN
  t = (a-0.5_dp) - 0.5_dp
  fn_val = gamln1(t)
  RETURN
END IF

IF (a < 10.0_dp) THEN
  n = a - 1.25_dp
  t = a
  w = 1.0_dp
  DO i = 1, n
    t = t - 1.0_dp
    w = t * w
  END DO
  fn_val = gamln1(t-1.0) + LOG(w)
  RETURN
END IF

t = (1.0/a) ** 2
w = (((((c5*t + c4)*t + c3)*t + c2)*t + c1)*t + c0) / a
fn_val = (d+w) + (a-0.5) * (LOG(a)-1.0)
RETURN
END FUNCTION gamln



FUNCTION algdiv(a, b) RESULT(fn_val)
!-----------------------------------------------------------------------

!     COMPUTATION OF LN(GAMMA(B)/GAMMA(A+B)) WHEN B >= 8

!                         --------

!     IN THIS ALGORITHM, DEL(X) IS THE FUNCTION DEFINED BY
!     LN(GAMMA(X)) = (X - 0.5)*LN(X) - X + 0.5*LN(2*PI) + DEL(X).

!-----------------------------------------------------------------------
IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(15, 100)

REAL (dp), INTENT(IN) :: a, b
REAL (dp)             :: fn_val

! Local variables
REAL (dp), PARAMETER  :: c0 = .833333333333333D-01, c1 = -.277777777760991D-02,  &
                    c2 = .793650666825390D-03, c3 = -.595202931351870D-03,  &
                    c4 = .837308034031215D-03, c5 = -.165322962780713D-02
REAL (dp) :: c, d, h, s3, s5, s7, s9, s11, t, u, v, w, x, x2
!------------------------
IF (a > b) THEN
  h = b / a
  c = 1.0 / (1.0_dp + h)
  x = h / (1.0_dp + h)
  d = a + (b - 0.5_dp)
ELSE
  h = a / b
  c = h / (1.0_dp + h)
  x = 1.0 / (1.0_dp + h)
  d = b + (a - 0.5_dp)
END IF

!                SET SN = (1 - X**N)/(1 - X)

x2 = x * x
s3 = 1.0 + (x + x2)
s5 = 1.0 + (x + x2*s3)
s7 = 1.0 + (x + x2*s5)
s9 = 1.0 + (x + x2*s7)
s11 = 1.0 + (x + x2*s9)

!                SET W = DEL(B) - DEL(A + B)

t = (1.0_dp/b) ** 2
w = ((((c5*s11*t + c4*s9)*t + c3*s7)*t + c2*s5)*t + c1*s3) * t + c0
w = w * (c/b)

!                    COMBINE THE RESULTS

u = d * alnrel(a/b)
v = a * (LOG(b) - 1.0_dp)
IF (u > v) THEN
  fn_val = (w-v) - u
  RETURN
END IF
fn_val = (w-u) - v
RETURN
END FUNCTION algdiv



FUNCTION gsumln(a, b) RESULT(fn_val)
!-----------------------------------------------------------------------
!          EVALUATION OF THE FUNCTION LN(GAMMA(A + B))
!          FOR 1 <= A <= 2  AND  1 <= B <= 2
!-----------------------------------------------------------------------
IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(15, 100)

REAL (dp), INTENT(IN) :: a, b
REAL (dp)             :: fn_val

! Local variables
REAL (dp) :: x

x = a + b - 2.0_dp
IF (x <= 0.25_dp) THEN
  fn_val = gamln1(1.0_dp + x)
  RETURN
END IF
IF (x <= 1.25_dp) THEN
  fn_val = gamln1(x) + alnrel(x)
  RETURN
END IF
fn_val = gamln1(x - 1.0_dp) + LOG(x*(1.0_dp + x))
RETURN
END FUNCTION gsumln



FUNCTION bcorr(a0, b0) RESULT(fn_val)
!-----------------------------------------------------------------------

!     EVALUATION OF  DEL(A0) + DEL(B0) - DEL(A0 + B0)  WHERE
!     LN(GAMMA(A)) = (A - 0.5)*LN(A) - A + 0.5*LN(2*PI) + DEL(A).
!     IT IS ASSUMED THAT A0 >= 8 AND B0 >= 8.

!-----------------------------------------------------------------------
IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(15, 100)

REAL (dp), INTENT(IN) :: a0, b0
REAL (dp)             :: fn_val

! Local variables
REAL (dp), PARAMETER :: c0 = .833333333333333D-01, c1 = -.277777777760991D-02,  &
                   c2 = .793650666825390D-03, c3 = -.595202931351870D-03,  &
                   c4 = .837308034031215D-03, c5 = -.165322962780713D-02
REAL (dp) :: a, b, c, h, s3, s5, s7, s9, s11, t, w, x, x2
!------------------------
a = MIN(a0,b0)
b = MAX(a0,b0)

h = a / b
c = h / (1.0_dp + h)
x = 1.0 / (1.0_dp + h)
x2 = x * x

!                SET SN = (1 - X**N)/(1 - X)

s3 = 1.0 + (x + x2)
s5 = 1.0 + (x + x2*s3)
s7 = 1.0 + (x + x2*s5)
s9 = 1.0 + (x + x2*s7)
s11 = 1.0 + (x + x2*s9)

!                SET W = DEL(B) - DEL(A + B)

t = (1.0_dp/b) ** 2
w = ((((c5*s11*t + c4*s9)*t + c3*s7)*t + c2*s5)*t + c1*s3) * t + c0
w = w * (c/b)

!                   COMPUTE  DEL(A) + W

t = (1.0_dp/a) ** 2
fn_val = (((((c5*t + c4)*t + c3)*t + c2)*t + c1)*t + c0) / a + w
RETURN
END FUNCTION bcorr



FUNCTION alnrel(a) RESULT(fn_val)
!-----------------------------------------------------------------------
!            EVALUATION OF THE FUNCTION LN(1 + A)
!-----------------------------------------------------------------------
IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(15, 100)

REAL (dp), INTENT(IN) :: a
REAL (dp)             :: fn_val

! Local variables
REAL (dp), PARAMETER  :: p1 = -.129418923021993D+01, p2 = .405303492862024D+00,  &
        p3 = -.178874546012214D-01, q1 = -.162752256355323D+01,   &
        q2 =  .747811014037616D+00, q3 = -.845104217945565D-01,  &
        zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp, two = 2.0_dp
REAL (dp) :: t, t2, w, x
!--------------------------
IF (ABS(a) <= 0.375_dp) THEN
  t = a / (a + two)
  t2 = t * t
  w = (((p3*t2 + p2)*t2 + p1)*t2 + one) / (((q3*t2 + q2)*t2 + q1)*t2 + one)
  fn_val = two * t * w
  RETURN
END IF

x = one + a
IF (a < zero) x = (a + half) + half
fn_val = LOG(x)
RETURN
END FUNCTION alnrel



FUNCTION gamln1(a) RESULT(fn_val)
!-----------------------------------------------------------------------
!     EVALUATION OF LN(GAMMA(1 + A)) FOR -0.2 <= A <= 1.25
!-----------------------------------------------------------------------
IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(15, 100)

REAL (dp), INTENT(IN) :: a
REAL (dp)             :: fn_val

! Local variables
REAL (dp), PARAMETER  :: p0 = .577215664901533D+00, p1 = .844203922187225D+00,  &
        p2 = -.168860593646662D+00, p3 = -.780427615533591D+00,   &
        p4 = -.402055799310489D+00, p5 = -.673562214325671D-01,   &
        p6 = -.271935708322958D-02, q1 = .288743195473681D+01,   &
        q2  = .312755088914843D+01, q3 = .156875193295039D+01,   &
        q4  = .361951990101499D+00, q5 = .325038868253937D-01,   &
        q6  = .667465618796164D-03, r0 = .422784335098467D+00,   &
        r1  = .848044614534529D+00, r2 = .565221050691933D+00,   &
        r3  = .156513060486551D+00, r4 = .170502484022650D-01,   &
        r5  = .497958207639485D-03, s1 = .124313399877507D+01,   &
        s2  = .548042109832463D+00, s3 = .101552187439830D+00,   &
        s4  = .713309612391000D-02, s5 = .116165475989616D-03
REAL (dp)  :: w, x
!----------------------
IF (a < 0.6_dp) THEN
  w = ((((((p6*a + p5)*a + p4)*a + p3)*a + p2)*a + p1)*a + p0) /   &
      ((((((q6*a + q5)*a + q4)*a + q3)*a + q2)*a + q1)*a + 1.0)
  fn_val = -a * w
  RETURN
END IF

x = (a - 0.5_dp) - 0.5_dp
w = (((((r5*x + r4)*x + r3)*x + r2)*x + r1)*x + r0) / (((((s5*x + s4)*x +  &
    s3)*x + s2)*x + s1)*x + 1.0_dp)
fn_val = x * w
RETURN
END FUNCTION gamln1

END MODULE LogBeta

MODULE Incomplete_Gamma
   use m_qtlmap_constant
!USE constants_NSWC
IMPLICIT NONE


CONTAINS


SUBROUTINE gratio(a, x, ans, qans, ind)
!-----------------------------------------------------------------------

!     EVALUATION OF THE INCOMPLETE GAMMA RATIO FUNCTIONS
!                   P(A,X) AND Q(A,X)

!                     ----------

!  IT IS ASSUMED THAT A AND X ARE NONNEGATIVE, WHERE A AND X ARE NOT BOTH 0.

!  ANS AND QANS ARE VARIABLES.  GRATIO ASSIGNS ANS THE VALUE P(A,X)
!  AND QANS THE VALUE Q(A,X). IND MAY BE ANY INTEGER.
!  IF IND = 0 THEN THE USER IS REQUESTING AS MUCH ACCURACY AS POSSIBLE
!  (UP TO 14 SIGNIFICANT DIGITS).  OTHERWISE, IF IND = 1 THEN ACCURACY
!  IND = 1 THEN ACCURACY IS REQUESTED TO WITHIN 1 UNIT OF THE 6-TH SIGNIFICANT
!  DIGIT, AND IF IND .NE. 0,1 THEN ACCURACY IS REQUESTED TO WITHIN 1 UNIT
!  OF THE 3RD SIGNIFICANT DIGIT.

!  ERROR RETURN ...

!     ANS IS ASSIGNED THE VALUE 2 WHEN A OR X IS NEGATIVE,
!  WHEN A*X = 0, OR WHEN P(A,X) AND Q(A,X) ARE INDETERMINANT.
!  P(A,X) AND Q(A,X) ARE COMPUTATIONALLY INDETERMINANT WHEN
!  X IS EXCEEDINGLY CLOSE TO A AND A IS EXTREMELY LARGE.

!--------------------------------------------------------------------
!  WRITTEN BY ALFRED H. MORRIS, JR.
!     NAVAL SURFACE WARFARE CENTER
!     DAHLGREN, VIRGINIA
!  REVISED ... DEC 1991
!-------------------------
REAL(dp), INTENT(IN)   :: a, x
INTEGER, INTENT(IN)    :: ind
REAL(dp), INTENT(OUT)  :: ans, qans

! Local variables
REAL(dp), PARAMETER :: acc0(3) = (/ 5.d-15, 5.d-7, 5.d-4 /),  &
        big(3) = (/ 25.0_dp, 14.0_dp , 10.0_dp /),  &
        e0(3) = (/ .25D-3, .25D-1, 0.14_dp /),  &
        x0(3) = (/ 31.0_dp, 17.0_dp, 9.7_dp /), ALOG10 = 2.30258509299405_dp,  &
        rt2pin = .398942280401433_dp, rtpi = 1.77245385090552_dp,  &
        d00 = -.333333333333333_dp, d10 = -.185185185185185D-02,   &
        d20 =  .413359788359788D-02, d30 = .649434156378601D-03,   &
        d40 = -.861888290916712D-03, d50 = -.336798553366358D-03,   &
        d60 =  .531307936463992D-03, d70 = .344367606892378D-03,   &
        d80 = -.652623918595309D-03
REAL(dp), PARAMETER :: a0(4) = (/ -.231272501940775D-02, -.335378520024220D-01,  &
                           -.159840143443990_dp, -.333333333333333_dp /),  &
        a1(4) = (/ -.398783924370770D-05, -.587926036018402D-03,  &
                   -.491687131726920D-02, -.185185185184291D-02 /),  &
        a2(2) = (/ .669564126155663D-03, .413359788442192D-02 /),  &
        a3(2) = (/ .810586158563431D-03, .649434157619770D-03 /),  &
        a4(2) = (/ -.105014537920131D-03, -.861888301199388D-03 /),  &
        a5(2) = (/ -.435211415445014D-03, -.336806989710598D-03 /),  &
        a6(2) = (/ -.182503596367782D-03, .531279816209452D-03 /),  &
        a7(2) = (/ .443219646726422D-03, .344430064306926D-03 /),  &
        a8(2) = (/ .878371203603888D-03, -.686013280418038D-03 /)
REAL(dp), PARAMETER :: b0(6) = (/ .633763414209504D-06, -.939001940478355D-05,  &
        .239521354917408D-02, .376245718289389D-01, .238549219145773_dp,  &
        .729520430331981_dp /),  &
        b1(4) = (/ .386325038602125D-02,  &
        .506042559238939D-01, .283344278023803_dp, .780110511677243_dp /),  &
        b2(5) = (/ -.421924263980656D-03, .650837693041777D-02,  &
        .682034997401259D-01, .339173452092224_dp, .810647620703045_dp /),  &
        b3(5) = (/ -.632276587352120D-03, .905375887385478D-02,  &
        .906610359762969D-01, .406288930253881_dp, .894800593794972_dp /),  &
        b4(4) = (/ .322609381345173D-01, .178295773562970_dp,  &
        .591353097931237_dp, .103151890792185D+01 /),  &
        b5(3) = (/ .178716720452422_dp, .600380376956324_dp,  &
        .108515217314415D+01 /),  &
        b6(2) = (/ .345608222411837_dp, .770341682526774_dp /),  &
        b7(2) = (/ .821824741357866_dp, .115029088777769D+01 /)
REAL(dp) :: d0(6) = (/ .833333333333333D-01, -.148148148148148D-01,  &
        .115740740740741D-02, .352733686067019D-03, -.178755144032922D-03,  &
        .391926317852244D-04 /), d1(4) = (/ -.347222222222222D-02,  &
        .264550264550265D-02, -.990226337448560D-03, .205761316872428D-03 /),  &
        d2(2) = (/ -.268132716049383D-02, .771604938271605D-03 /),  &
        d3(2) = (/ .229472093621399D-03, -.469189494395256D-03 /),  &
        d4(1) = (/ .784039221720067D-03 /), d5(1) = (/ -.697281375836586D-04 /),  &
        d6(1) = (/ -.592166437353694D-03 /)
REAL(dp) :: acc, amn, apn, a2n, a2nm1, b2n, b2nm1, c, c0, c1, c2, c3, c4, c5, &
            c6, c7, c8, e, g, h, j, l, r, rta, rtx, s, sum, t, tol, twoa, u,  &
            w, wk(20), y, z
INTEGER  :: i, iop, m, n, nl1
!-------------------------

!     ****** E IS A MACHINE DEPENDENT CONSTANT. E IS THE SMALLEST
!            FLOATING POINT NUMBER FOR WHICH 1.0 + E > 1.0 .

e = EPSILON(1.0_dp)

!-------------------------
IF (a < 0.0 .OR. x < 0.0) GO TO 320
IF (a == 0.0 .AND. x == 0.0) GO TO 320
IF (a*x == 0.0) GO TO 310

iop = ind + 1
IF (iop /= 1 .AND. iop /= 2) iop = 3
acc = MAX(acc0(iop),e)

!            SELECT THE APPROPRIATE ALGORITHM

IF (a < 1.0) THEN
  IF (a == 0.5) GO TO 290
  IF (x < 1.1) GO TO 100
  r = drcomp(a,x)
  IF (r == 0.0) GO TO 280
  GO TO 160
END IF

IF (a < big(iop)) THEN
  IF (a > x .OR. x >= x0(iop)) GO TO 10
  twoa = a + a
  m = INT(twoa)
  IF (twoa /= REAL(m)) GO TO 10
  i = m / 2
  IF (a == REAL(i)) GO TO 130
  GO TO 140
END IF

l = x / a
IF (l == 0.0) GO TO 270
s = 0.5 + (0.5-l)
z = drlog(l)
IF (z >= 700.0/a) GO TO 300
y = a * z
rta = SQRT(a)
IF (ABS(s) <= e0(iop)/rta) GO TO 230
IF (ABS(s) <= 0.4) GO TO 180

10 r = drcomp(a,x)
IF (r == 0.0) GO TO 310
IF (x > MAX(a, ALOG10)) THEN
  IF (x < x0(iop)) GO TO 160
ELSE

!                 TAYLOR SERIES FOR P/R

  apn = a + 1.0
  t = x / apn
  wk(1) = t
  DO n = 2, 20
    apn = apn + 1.0
    t = t * (x/apn)
    IF (t <= 1.D-3) GO TO 30
    wk(n) = t
  END DO
  n = 20

  30 sum = t
  tol = 0.5 * acc
  40 apn = apn + 1.0
  t = t * (x/apn)
  sum = sum + t
  IF (t > tol) GO TO 40

  nl1 = n - 1
  DO m = 1, nl1
    n = n - 1
    sum = sum + wk(n)
  END DO
  ans = (r/a) * (1.0+sum)
  qans = 0.5 + (0.5-ans)
  RETURN
END IF

!                 ASYMPTOTIC EXPANSION

amn = a - 1.0
t = amn / x
wk(1) = t
DO n = 2, 20
  amn = amn - 1.0
  t = t * (amn/x)
  IF (ABS(t) <= 1.D-3) GO TO 70
  wk(n) = t
END DO
n = 20

70 sum = t
80 IF (ABS(t) >= acc) THEN
  amn = amn - 1.0
  t = t * (amn/x)
  sum = sum + t
  GO TO 80
END IF

nl1 = n - 1
DO m = 1, nl1
  n = n - 1
  sum = sum + wk(n)
END DO
qans = (r/x) * (1.0+sum)
ans = 0.5 + (0.5-qans)
RETURN

!             TAYLOR SERIES FOR P(A,X)/X**A

100 l = 3.0
c = x
sum = x / (a+3.0)
tol = 3.0 * acc / (a+1.0)
110 l = l + 1.0
c = -c * (x/l)
t = c / (a+l)
sum = sum + t
IF (ABS(t) > tol) GO TO 110
j = a * x * ((sum/6.0 - 0.5/(a+2.0))*x + 1.0/(a+1.0))

z = a * LOG(x)
h = dgam1(a)
g = 1.0 + h
IF (x >= 0.25) THEN
  IF (a < x/2.59) GO TO 120
ELSE
  IF (z > -.13394) GO TO 120
END IF

w = EXP(z)
ans = w * g * (0.5+(0.5-j))
qans = 0.5 + (0.5-ans)
RETURN

120 l = drexp(z)
w = 0.5 + (0.5+l)
qans = (w*j-l) * g - h
IF (qans < 0.0) GO TO 280
ans = 0.5 + (0.5-qans)
RETURN

!             FINITE SUMS FOR Q WHEN A >= 1
!                 AND 2*A IS AN INTEGER

130 sum = EXP(-x)
t = sum
n = 1
c = 0.0
GO TO 150

140 rtx = SQRT(x)
sum = derfc1(0,rtx)
t = EXP(-x) / (rtpi*rtx)
n = 0
c = -0.5

150 IF (n /= i) THEN
  n = n + 1
  c = c + 1.0
  t = (x*t) / c
  sum = sum + t
  GO TO 150
END IF
qans = sum
ans = 0.5 + (0.5-qans)
RETURN

!              CONTINUED FRACTION EXPANSION

160 tol = MAX(8.0*e,4.0*acc)
a2nm1 = 1.0
a2n = 1.0
b2nm1 = x
b2n = x + (1.0-a)
c = 1.0
170 a2nm1 = x * a2n + c * a2nm1
b2nm1 = x * b2n + c * b2nm1
c = c + 1.0
t = c - a
a2n = a2nm1 + t * a2n
b2n = b2nm1 + t * b2n

a2nm1 = a2nm1 / b2n
b2nm1 = b2nm1 / b2n
a2n = a2n / b2n
b2n = 1.0
IF (ABS(a2n-a2nm1/b2nm1) >= tol*a2n) GO TO 170

qans = r * a2n
ans = 0.5 + (0.5-qans)
RETURN

180 IF (ABS(s) <= 2.0*e .AND. a*e*e > 3.28D-3) GO TO 320
c = EXP(-y)
w = 0.5 * derfc1(1,SQRT(y))
u = 1.0 / a
z = SQRT(z+z)
IF (l < 1.0) z = -z
IF (iop == 2) THEN
  GO TO 200
ELSE IF (iop > 2) THEN
  GO TO 210
END IF

IF (ABS(s) <= 1.D-3) GO TO 240

!            USING THE MINIMAX APPROXIMATIONS

c0 = (((a0(1)*z + a0(2))*z + a0(3))*z + a0(4)) / ((((((b0(1)*z + b0(2))*z +  &
     b0(3))*z + b0(4))*z + b0(5))*z + b0(6))*z + 1.0)
c1 = (((a1(1)*z + a1(2))*z + a1(3))*z + a1(4)) / ((((b1(1)*z + b1(2))*z +   &
     b1(3))*z + b1(4))*z + 1.0)
c2 = (a2(1)*z + a2(2)) / (((((b2(1)*z + b2(2))*z + b2(3))*z + b2(4))*z +  &
     b2(5))*z + 1.0)
c3 = (a3(1)*z + a3(2)) / (((((b3(1)*z + b3(2))*z + b3(3))*z + b3(4))*z +  &
     b3(5))*z + 1.0)
c4 = (a4(1)*z + a4(2)) / ((((b4(1)*z + b4(2))*z + b4(3))*z + b4(4))*z + 1.0)
c5 = (a5(1)*z + a5(2)) / (((b5(1)*z + b5(2))*z + b5(3))*z + 1.0)
c6 = (a6(1)*z + a6(2)) / ((b6(1)*z + b6(2))*z + 1.0)
c7 = (a7(1)*z + a7(2)) / ((b7(1)*z + b7(2))*z + 1.0)
c8 = a8(1) * z + a8(2)
t = (((((((c8*u + c7)*u + c6)*u + c5)*u + c4)*u + c3)*u + c2)*u + c1)*u + c0
GO TO 220

!                    TEMME EXPANSION

200 c0 = (((((d0(6)*z+d0(5))*z+d0(4))*z+d0(3))*z+d0(2))*z+d0(1)) * z +d00
c1 = (((d1(4)*z+d1(3))*z+d1(2))*z+d1(1)) * z + d10
c2 = d2(1) * z + d20
t = (c2*u+c1) * u + c0
GO TO 220

210 t = ((d0(3)*z+d0(2))*z+d0(1)) * z + d00

220 IF (l >= 1.0) THEN
  qans = c * (w+rt2pin*t/rta)
  ans = 0.5 + (0.5-qans)
  RETURN
END IF
ans = c * (w-rt2pin*t/rta)
qans = 0.5 + (0.5-ans)
RETURN

!               TEMME EXPANSION FOR L = 1

230 IF (a*e*e > 3.28D-3) GO TO 320
c = 0.5 + (0.5-y)
w = (0.5 - SQRT(y)*(0.5 + (0.5 - y/3.0))/rtpi) / c
u = 1.0 / a
z = SQRT(z+z)
IF (l < 1.0) z = -z
IF (iop == 2) THEN
  GO TO 250
ELSE IF (iop > 2) THEN
  GO TO 260
END IF

240 c0 = ((d0(3)*z+d0(2))*z+d0(1)) * z + d00
c1 = ((d1(3)*z+d1(2))*z+d1(1)) * z + d10
c2 = (d2(2)*z+d2(1)) * z + d20
c3 = (d3(2)*z+d3(1)) * z + d30
c4 = d4(1) * z + d40
c5 = d5(1) * z + d50
c6 = d6(1) * z + d60
t = (((((((d80*u+d70)*u+c6)*u+c5)*u+c4)*u+c3)*u+c2)*u+c1) * u + c0
GO TO 220

250 c0 = (d0(2)*z+d0(1)) * z + d00
c1 = d1(1) * z + d10
t = (d20*u+c1) * u + c0
GO TO 220

260 t = d0(1) * z + d00
GO TO 220

!                     SPECIAL CASES

270 ans = 0.0
qans = 1.0
RETURN

280 ans = 1.0
qans = 0.0
RETURN

290 IF (x < 0.25) THEN
  ans = derf(SQRT(x))
  qans = 0.5 + (0.5-ans)
  RETURN
END IF
qans = derfc1(0,SQRT(x))
ans = 0.5 + (0.5-qans)
RETURN

300 IF (ABS(s) <= 2.0*e) GO TO 320
310 IF (x <= a) GO TO 270
GO TO 280

!                     ERROR RETURN

320 ans = 2.0
RETURN
END SUBROUTINE gratio



SUBROUTINE gaminv(a, x, x0, p, q, ierr)
!--------------------------------------------------------------------

!          INVERSE INCOMPLETE GAMMA RATIO FUNCTION

!  GIVEN POSITIVE A, AND NONEGATIVE P AND Q WHERE P + Q = 1.
!  THEN X IS COMPUTED WHERE P(A,X) = P AND Q(A,X) = Q. SCHRODER
!  ITERATION IS EMPLOYED. THE ROUTINE ATTEMPTS TO COMPUTE X
!  TO 10 SIGNIFICANT DIGITS IF THIS IS POSSIBLE FOR THE
!  PARTICULAR COMPUTER ARITHMETIC BEING USED.

!                     ------------

!  X IS A VARIABLE. IF P = 0 THEN X IS ASSIGNED THE VALUE 0,
!  AND IF Q = 0 THEN X IS SET TO THE LARGEST FLOATING POINT
!  NUMBER AVAILABLE. OTHERWISE, GAMINV ATTEMPTS TO OBTAIN
!  A SOLUTION FOR P(A,X) = P AND Q(A,X) = Q. IF THE ROUTINE
!  IS SUCCESSFUL THEN THE SOLUTION IS STORED IN X.

!  X0 IS AN OPTIONAL INITIAL APPROXIMATION FOR X. IF THE USER DOES NOT
!  WISH TO SUPPLY AN INITIAL APPROXIMATION, THEN SET X0 <= 0.

!  IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.
!  WHEN THE ROUTINE TERMINATES, IERR HAS ONE OF THE FOLLOWING
!  VALUES ...

!    IERR =  0    THE SOLUTION WAS OBTAINED.  ITERATION WAS NOT USED.
!    IERR >  0    THE SOLUTION WAS OBTAINED.  IERR ITERATIONS WERE PERFORMED.
!    IERR = -2    (INPUT ERROR) A <= 0
!    IERR = -3    NO SOLUTION WAS OBTAINED. THE RATIO Q/A IS TOO LARGE.
!    IERR = -4    (INPUT ERROR) P OR Q IS NEGATIVE, OR P + Q .NE. 1.
!    IERR = -6    20 ITERATIONS WERE PERFORMED. THE MOST
!                 RECENT VALUE OBTAINED FOR X IS GIVEN.
!                 THIS CANNOT OCCUR IF X0 <= 0.
!    IERR = -7    ITERATION FAILED. NO VALUE IS GIVEN FOR X.
!                 THIS MAY OCCUR WHEN X IS APPROXIMATELY 0.
!    IERR = -8    A VALUE FOR X HAS BEEN OBTAINED, BUT THE ROUTINE IS NOT
!                 CERTAIN OF ITS ACCURACY.
!                 ITERATION CANNOT BE PERFORMED IN THIS CASE.
!                 IF X0 <= 0, THIS CAN OCCUR ONLY WHEN P OR Q IS
!                 APPROXIMATELY 0.  IF X0 IS POSITIVE THEN THIS CAN OCCUR
!                 WHEN A IS EXCEEDINGLY CLOSE TO X AND A IS EXTREMELY
!                 LARGE (SAY A >= 1.E20).

!--------------------------------------------------------------------
!  WRITTEN BY ALFRED H. MORRIS, JR.
!     NAVAL SURFACE WARFARE CENTER
!     DAHLGREN, VIRGINIA
!  REVISED ... JANUARY 1992
!------------------------
REAL(dp), INTENT(IN)   :: a, x0, p, q
REAL(dp), INTENT(OUT)  :: x
INTEGER, INTENT(OUT)   :: ierr

! Local variables
REAL(dp), PARAMETER :: ln10 = 2.302585_dp, bmin(2) = (/ 1.D-28, 1.D-13 /),  &
        emin(2) = (/ 2.D-03, 6.D-03 /), c = .577215664901533_dp, tol = 1.D-5
REAL(dp)  :: amax, amin, am1, ap1, ap2, ap3, apn, b, c1, c2, c3, c4, c5, d,  &
             e, eps, e2, g, h, pn, qg, qn, r, rta, s, sum, s2, t, u, w, xmin, &
             xn, y, z,zero,two
INTEGER   :: ier, iop

zero= 0.0
two=2.0
!------------------------
!     LN10 = LN(10)
!     C = EULER CONSTANT
!------------------------

!     ****** E AND XMIN ARE MACHINE DEPENDENT CONSTANTS. E IS THE
!            SMALLEST NUMBER FOR WHICH 1.0 + E > 1.0, AND XMIN
!            IS THE SMALLEST POSITIVE NUMBER.

e = EPSILON(1.0_dp)
xmin = TINY(1.0_dp)

!------------------------
x = 0.0
IF (a > 0.0) THEN
  IF (p < 0.0 .OR. q < 0.0) GO TO 120
  t = ((p+q)-0.5) - 0.5
  IF (ABS(t) > 5.0*MAX(e,1.D-15)) GO TO 120

  ierr = 0
  xmin = xmin / e
  IF ((p/e) > xmin) THEN
    IF ((q/e) <= xmin) GO TO 160
    IF (a == 1.0) GO TO 100

    e2 = e + e
    amax = 0.4D-10 / (e*e)
    eps = MAX(100.0*e,1.D-10)
    iop = 1
    IF (e > 1.D-10) iop = 2
    xn = x0
    IF (x0 <= 0.0) THEN

!        SELECTION OF THE INITIAL APPROXIMATION XN OF X
!                       WHEN A < 1

      IF (a > 1.0) GO TO 30
      g = dgamma(a+1.0)
      qg = q * g
      IF (qg == 0.0) GO TO 160
      b = qg / a
      IF (qg > 0.6*a) GO TO 20
      IF (a < 0.30 .AND. b >= 0.35) THEN
        t = EXP(-(b+c))
        u = t * EXP(t)
        xn = t * EXP(u)
        GO TO 50
      END IF

      IF (b >= 0.45) GO TO 20
      IF (b == 0.0) GO TO 160
      y = -LOG(b)
      s = 0.5 + (0.5-a)
      z = LOG(y)
      t = y - s * z
      IF (b >= 0.15) THEN
        xn = y - s * LOG(t) - LOG(1.0 + s/(t+1.0))
        GO TO 80
      END IF
      IF (b > 1.D-2) THEN
        u = ((t + 2.0*(3.0-a))*t + (2.0-a)*(3.0-a)) / ((t + (5.0-a))*t + 2.0)
        xn = y - s * LOG(t) - LOG(u)
        GO TO 80
      END IF
      10 c1 = -s * z
      c2 = -s * (1.0+c1)
      c3 = s * ((0.5*c1+(2.0-a))*c1+(2.5-1.5*a))
      c4 = -s * (((c1/3.0+(2.5-1.5*a))*c1+((a-6.0)*a+7.0))*c1 +  &
                ((11.0*a-46.0)*a+47.0)/6.0)
      c5 = -s * ((((-c1/4.0+(11.0*a-17.0)/6.0)*c1+((-3.0*a+13.0)*a  &
           -13.0))*c1+0.5*(((2.0*a-25.0)*a+72.0)*a-61.0))*c1 +  &
           (((25.0*a-195.0)*a+477.0)*a-379.0)/12.0)
      xn = ((((c5/y+c4)/y+c3)/y+c2)/y+c1) + y
      IF (a > 1.0) GO TO 80
      IF (b > bmin(iop)) GO TO 80
      x = xn
      RETURN

      20 IF (b*q <= 1.D-8) THEN
        xn = EXP(-(q/a+c))
      ELSE
        IF (p > 0.9) THEN
          xn = EXP((dlnrel(-q) + dgmln1(a))/a)
        ELSE
          xn = EXP(LOG(p*g)/a)
        END IF
      END IF

      IF (xn == 0.0) GO TO 110
      t = 0.5 + (0.5 - xn/(a+1.0))
      xn = xn / t
      GO TO 50

!        SELECTION OF THE INITIAL APPROXIMATION XN OF X
!                       WHEN A > 1

      30 t = p - 0.5
      IF (q < 0.5) t = 0.5 - q
      CALL dpni(p, q, t, s, ier)
      IF (ier /= 0) WRITE(*, *) '** Error in call to PNI from GAMINV **'

      rta = SQRT(a)
      s2 = s * s
      xn = (((12.0*s2 - 243.0)*s2 - 923.0)*s2 + 1472.0) / 204120.0
      xn = (xn/a + s*((9.0*s2 + 256.0)*s2 - 433.0)/(38880.0*rta)) -   &
           ((3.0*s2 + 7.0)*s2 - 16.0) / 810.0
      xn = a + s * rta + (s2-1.0) / 3.0 + s * (s2-7.0) / (36.0*rta) + xn / a
      xn = MAX(xn, zero)

      amin = 20.0
      IF (e < 1.D-8) amin = 250.0
      IF (a >= amin) THEN
        x = xn
        d = 0.5 + (0.5-x/a)
        IF (ABS(d) <= 1.D-1) RETURN
      END IF

      IF (p > 0.5) THEN
        IF (xn < 3.0*a) GO TO 80
        w = LOG(q)
        y = -(w + dgamln(a))
        d = MAX(two,a*(a-1.0))
        IF (y >= ln10*d) THEN
          s = 1.0 - a
          z = LOG(y)
          GO TO 10
        END IF
        t = a - 1.0
        xn = y + t * LOG(xn) - dlnrel(-t/(xn+1.0))
        xn = y + t * LOG(xn) - dlnrel(-t/(xn+1.0))
        GO TO 80
      END IF

      ap1 = a + 1.0
      IF (xn > 0.70*ap1) GO TO 60
      w = LOG(p) + dgamln(ap1)
      IF (xn <= 0.15*ap1) THEN
        ap2 = a + 2.0
        ap3 = a + 3.0
        x = EXP((w+x)/a)
        x = EXP((w+x-LOG(1.0+(x/ap1)*(1.0+x/ap2)))/a)
        x = EXP((w+x-LOG(1.0+(x/ap1)*(1.0+x/ap2)))/a)
        x = EXP((w+x-LOG(1.0+(x/ap1)*(1.0+(x/ap2)*(1.0+x/ap3))))/a)
        xn = x
        IF (xn <= 1.D-2*ap1) THEN
          IF (xn <= emin(iop)*ap1) RETURN
          GO TO 60
        END IF
      END IF

      apn = ap1
      t = xn / apn
      sum = 1.0 + t
      40 apn = apn + 1.0
      t = t * (xn/apn)
      sum = sum + t
      IF (t > 1.D-4) GO TO 40
      t = w - LOG(sum)
      xn = EXP((xn+t)/a)
      xn = xn * (1.0-(a*LOG(xn)-xn-t)/(a-xn))
      GO TO 60
    END IF

!                 SCHRODER ITERATION USING P

    50 IF (p > 0.5) GO TO 80
    60 IF (p <= xmin) GO TO 150
    am1 = (a-0.5) - 0.5
    70 IF (a > amax) THEN
      d = 0.5 + (0.5-xn/a)
      IF (ABS(d) <= e2) GO TO 150
    END IF

    IF (ierr >= 20) GO TO 130
    ierr = ierr + 1
    CALL gratio(a,xn,pn,qn,0)
    IF (pn == 0.0.OR.qn == 0.0) GO TO 150
    r = drcomp(a,xn)
    IF (r < xmin) GO TO 150
    t = (pn-p) / r
    w = 0.5 * (am1-xn)
    IF (ABS(t) > 0.1.OR.ABS(w*t) > 0.1) THEN
      x = xn * (1.0-t)
      IF (x <= 0.0) GO TO 140
      d = ABS(t)
    ELSE

      h = t * (1.0+w*t)
      x = xn * (1.0-h)
      IF (x <= 0.0) GO TO 140
      IF (ABS(w) >= 1.0 .AND. ABS(w)*t*t <= eps) RETURN
      d = ABS(h)
    END IF
    xn = x
    IF (d > tol) GO TO 70
    IF (d <= eps) RETURN
    IF (ABS(p-pn) <= tol*p) RETURN
    GO TO 70

!                 SCHRODER ITERATION USING Q

    80 IF (q <= xmin) GO TO 150
    am1 = (a-0.5) - 0.5
    90 IF (a > amax) THEN
      d = 0.5 + (0.5-xn/a)
      IF (ABS(d) <= e2) GO TO 150
    END IF

    IF (ierr >= 20) GO TO 130
    ierr = ierr + 1
    CALL gratio(a,xn,pn,qn,0)
    IF (pn == 0.0 .OR. qn == 0.0) GO TO 150
    r = drcomp(a,xn)
    IF (r < xmin) GO TO 150
    t = (q-qn) / r
    w = 0.5 * (am1-xn)
    IF (ABS(t) > 0.1 .OR. ABS(w*t) > 0.1) THEN
      x = xn * (1.0-t)
      IF (x <= 0.0) GO TO 140
      d = ABS(t)
    ELSE

      h = t * (1.0+w*t)
      x = xn * (1.0-h)
      IF (x <= 0.0) GO TO 140
      IF (ABS(w) >= 1.0 .AND. ABS(w)*t*t <= eps) RETURN
      d = ABS(h)
    END IF
    xn = x
    IF (d > tol) GO TO 90
    IF (d <= eps) RETURN
    IF (ABS(q-qn) <= tol*q) RETURN
    GO TO 90
  END IF

!                       SPECIAL CASES

  ierr = -8
  RETURN

  100 IF (q >= 0.9) THEN
    x = -dlnrel(-p)
    RETURN
  END IF
  x = -LOG(q)
  RETURN
END IF

!                       ERROR RETURN

ierr = -2
RETURN

110 ierr = -3
RETURN

120 ierr = -4
RETURN

130 ierr = -6
RETURN

140 ierr = -7
RETURN

150 x = xn
ierr = -8
RETURN

160 x = HUGE(1.0)
ierr = -8
RETURN
END SUBROUTINE gaminv



FUNCTION derf(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!        REAL (dp) EVALUATION OF THE ERROR FUNCTION
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN)  :: x
REAL (dp)              :: fn_val

! Local variables
REAL (dp)  :: ax, t, w
INTEGER    :: i, k
REAL (dp), PARAMETER :: a(21) = (/ .1283791670955125738961589031215_dp,  &
        -.3761263890318375246320529677070_dp,  &
        .1128379167095512573896158902931_dp,  &
        -.2686617064513125175943235372542D-01,  &
        .5223977625442187842111812447877D-02,  &
        -.8548327023450852832540164081187D-03,  &
        .1205533298178966425020717182498D-03,  &
        -.1492565035840625090430728526820D-04,  &
        .1646211436588924261080723578109D-05,  &
        -.1636584469123468757408968429674D-06,  &
        .1480719281587021715400818627811D-07,  &
        -.1229055530145120140800510155331D-08,  &
        .9422759058437197017313055084212D-10,  &
        -.6711366740969385085896257227159D-11,  &
        .4463222608295664017461758843550D-12,  &
        -.2783497395542995487275065856998D-13,  &
        .1634095572365337143933023780777D-14,  &
        -.9052845786901123985710019387938D-16,  &
        .4708274559689744439341671426731D-17,  &
        -.2187159356685015949749948252160D-18,  &
        .7043407712019701609635599701333D-20 /)
!-------------------------------

!                     ABS(X) <= 1

ax = ABS(x)
IF (ax <= 1._dp) THEN
  t = x * x
  w = a(21)
  DO i = 1, 20
    k = 21 - i
    w = t * w + a(k)
  END DO
  fn_val = x * (1._dp+w)
  RETURN
END IF

!                     ABS(X) > 1

IF (ax < 8.5_dp) THEN
  fn_val = 0.5_dp + (0.5_dp - EXP(-x*x)*derfc0(ax))
  IF (x < 0._dp) fn_val = -fn_val
  RETURN
END IF

!                 LIMIT VALUE FOR LARGE X

fn_val = SIGN(1._dp,x)
RETURN
END FUNCTION derf



FUNCTION derfc1(ind, x) RESULT(fn_val)
!--------------------------------------------------------------------

!      EVALUATION OF THE COMPLEMENTARY ERROR FUNCTION

!       DERFC1(IND,X) = ERFC(X)           IF IND = 0
!       DERFC1(IND,X) = EXP(X*X)*ERFC(X)  OTHERWISE

!--------------------------------------------------------------------
INTEGER, INTENT(IN)    :: ind
REAL (dp), INTENT(IN)  :: x
REAL (dp)              :: fn_val

! Local variables
REAL (dp)  :: ax, t, w
INTEGER    :: i, k
REAL (dp), PARAMETER :: a(21) = (/ .1283791670955125738961589031215_dp,  &
        -.3761263890318375246320529677070_dp,  &
        .1128379167095512573896158902931_dp,  &
        -.2686617064513125175943235372542D-01,  &
        .5223977625442187842111812447877D-02,  &
        -.8548327023450852832540164081187D-03,  &
        .1205533298178966425020717182498D-03,  &
        -.1492565035840625090430728526820D-04,  &
        .1646211436588924261080723578109D-05,  &
        -.1636584469123468757408968429674D-06,  &
        .1480719281587021715400818627811D-07,  &
        -.1229055530145120140800510155331D-08,  &
        .9422759058437197017313055084212D-10,  &
        -.6711366740969385085896257227159D-11,  &
        .4463222608295664017461758843550D-12,  &
        -.2783497395542995487275065856998D-13,  &
        .1634095572365337143933023780777D-14,  &
        -.9052845786901123985710019387938D-16,  &
        .4708274559689744439341671426731D-17,  &
        -.2187159356685015949749948252160D-18,  &
        .7043407712019701609635599701333D-20 /)
!-------------------------------

!                     ABS(X) <= 1

ax = ABS(x)
IF (ax <= 1._dp) THEN
  t = x * x
  w = a(21)
  DO i = 1, 20
    k = 21 - i
    w = t * w + a(k)
  END DO
  fn_val = 0.5_dp + (0.5_dp-x*(1._dp+w))
  IF (ind /= 0) fn_val = EXP(t) * fn_val
  RETURN
END IF

!                       X < -1

IF (x <= 0._dp) THEN
  IF (x < -8.3_dp) GO TO 20
  IF (ind /= 0) THEN
    fn_val = 2._dp * EXP(x*x) - derfc0(ax)
    RETURN
  END IF
  fn_val = 2._dp - EXP(-x*x) * derfc0(ax)
  RETURN
END IF

!                       X > 1

IF (ind /= 0) THEN
  fn_val = derfc0(x)
  RETURN
END IF
fn_val = 0._dp
IF (x > 100._dp) RETURN
t = x * x
IF (t > -dxparg(1)) RETURN
fn_val = EXP(-t) * derfc0(x)
RETURN

!             LIMIT VALUE FOR LARGE NEGATIVE X

20 fn_val = 2._dp
IF (ind /= 0) fn_val = 2._dp * EXP(x*x)
RETURN
END FUNCTION derfc1


FUNCTION derfc0(x) RESULT(fn_val)
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN)  :: x
REAL (dp)              :: fn_val

!        EVALUATION OF EXP(X**2)*ERFC(X) FOR X >= 1

!--------------------------------------------------------------------
!  WRITTEN BY ALFRED H. MORRIS, JR.
!     NAVAL SURFACE WARFARE CENTER
!     DAHLGREN, VIRGINIA
!     APRIL 1992
!-------------------------------
REAL (dp)            :: t, u, v, z
REAL (dp), PARAMETER :: rpinv = .56418958354775628694807945156077259_dp
REAL (dp), PARAMETER :: p0 = .16506148041280876191828601D-03,  &
                        p1 =  .15471455377139313353998665D-03,  &
                        p2 =  .44852548090298868465196794D-04,  &
                        p3 = -.49177280017226285450486205D-05,  &
                        p4 = -.69353602078656412367801676D-05,  &
                        p5 = -.20508667787746282746857743D-05,  &
                        p6 = -.28982842617824971177267380D-06,  &
                        p7 = -.17272433544836633301127174D-07,  &
                        q1 =  .16272656776533322859856317D+01,  &
                        q2 =  .12040996037066026106794322D+01,  &
                        q3 =  .52400246352158386907601472_dp,  &
                        q4 =  .14497345252798672362384241_dp,  &
                        q5 =  .25592517111042546492590736D-01,  &
                        q6 =  .26869088293991371028123158D-02,  &
                        q7 =  .13133767840925681614496481D-03
REAL (dp), PARAMETER :: r0 =  .145589721275038539045668824025_dp,  &
                        r1 = -.273421931495426482902320421863_dp,  &
                        r2 =  .226008066916621506788789064272_dp,  &
                        r3 = -.163571895523923805648814425592_dp,  &
                        r4 =  .102604312032193978662297299832_dp,  &
                        r5 = -.548023266949835519254211506880D-01,  &
                        r6 =  .241432239725390106956523668160D-01,  &
                        r7 = -.822062115403915116036874169600D-02,  &
                        r8 =  .180296241564687154310619200000D-02
REAL (dp), PARAMETER :: a0 = -.45894433406309678202825375D-03,   &
                        a1 = -.12281298722544724287816236D-01,  &
                        a2 = -.91144359512342900801764781D-01,  &
                        a3 = -.28412489223839285652511367D-01,  &
                        a4 =  .14083827189977123530129812D+01,  &
                        a5 =  .11532175281537044570477189D+01,  &
                        a6 = -.72170903389442152112483632D+01,  &
                        a7 = -.19685597805218214001309225D+01,  &
                        a8 =  .93846891504541841150916038D+01,  &
                        b1 =  .25136329960926527692263725D+02,  &
                        b2 =  .15349442087145759184067981D+03,  &
                        b3 = -.29971215958498680905476402D+03,  &
                        b4 = -.33876477506888115226730368D+04,  &
                        b5 =  .28301829314924804988873701D+04,  &
                        b6 =  .22979620942196507068034887D+05,  &
                        b7 = -.24280681522998071562462041D+05,  &
                        b8 = -.36680620673264731899504580D+05,  &
                        b9 =  .42278731622295627627042436D+05,  &
                        b10=  .28834257644413614344549790D+03,  &
                        b11=  .70226293775648358646587341D+03
REAL (dp), PARAMETER :: c0 = -.7040906288250128001000086D-04,   &
                        c1 = -.3858822461760510359506941D-02,  &
                        c2 = -.7708202127512212359395078D-01,  &
                        c3 = -.6713655014557429480440263_dp,  &
                        c4 = -.2081992124162995545731882D+01,  &
                        c5 =  .2898831421475282558867888D+01,  &
                        c6 =  .2199509380600429331650192D+02,  &
                        c7 =  .2907064664404115316722996D+01,  &
                        c8 = -.4766208741588182425380950D+02,  &
                        d1 =  .5238852785508439144747174D+02,  &
                        d2 =  .9646843357714742409535148D+03,  &
                        d3 =  .7007152775135939601804416D+04,  &
                        d4 =  .8515386792259821780601162D+04,  &
                        d5 = -.1002360095177164564992134D+06,  &
                        d6 = -.2065250031331232815791912D+06,  &
                        d7 =  .5695324805290370358175984D+06,  &
                        d8 =  .6589752493461331195697873D+06,  &
                        d9 = -.1192930193156561957631462D+07
REAL (dp), PARAMETER :: e0 = .540464821348814822409610122136_dp,  &
                        e1 = -.261515522487415653487049835220D-01, &
                        e2 = -.288573438386338758794591212600D-02, &
                        e3 = -.529353396945788057720258856000D-03
REAL (dp), PARAMETER :: s1 = .75000000000000000000_dp,   &
        s2  = -.18750000000000000000D+01, s3  = .65625000000000000000D+01,  &
        s4  = -.29531250000000000000D+02, s5  = .16242187500000000000D+03,  &
        s6  = -.10557421875000000000D+04, s7  = .79180664062500000000D+04,  &
        s8  = -.67303564453125000000D+05, s9  = .63938386230468750000D+06,  &
        s10 = -.67135305541992187500D+07, s11 = .77205601373291015625D+08
!-------------------------------
!     RPINV = 1/SQRT(PI)
!-------------------------------

!                     1 <= X <= 2

IF (x <= 2._dp) THEN
  u = ((((((p7*x + p6)*x + p5)*x + p4)*x + p3)*x + p2)*x + p1) * x + p0
  v = ((((((q7*x + q6)*x + q5)*x + q4)*x + q3)*x + q2)*x + q1) * x + 1._dp
  t = (x-3.75_dp) / (x+3.75_dp)
  fn_val = (((((((((u/v)*t + r8)*t + r7)*t + r6)*t + r5)*t + r4)*t + r3)*t + &
           r2)*t + r1) * t + r0
  RETURN
END IF

!                     2 < X <= 4

IF (x <= 4._dp) THEN
  z = 1._dp / (2.5_dp + x*x)
  u = (((((((a8*z + a7)*z + a6)*z + a5)*z + a4)*z + a3)*z + a2)*z + a1) * z + a0
  v = ((((((((((b11*z + b10)*z + b9)*z + b8)*z + b7)*z + b6)*z + b5)*z +  &
      b4)*z + b3)*z + b2)*z + b1) * z + 1._dp
  t = 13._dp * z - 1._dp
  fn_val = ((((u/v)*t + e2)*t + e1)*t + e0) / x
  RETURN
END IF

!                     4 < X < 50

IF (x < 50._dp) THEN
  z = 1._dp / (2.5_dp + x*x)
  u = (((((((c8*z + c7)*z + c6)*z + c5)*z + c4)*z + c3)*z + c2)*z + c1) * z + &
      c0
  v = ((((((((d9*z + d8)*z + d7)*z + d6)*z + d5)*z + d4)*z + d3)*z + d2)*z +  &
      d1)*z + 1._dp
  t = 13._dp * z - 1._dp
  fn_val = (((((u/v)*t + e3)*t + e2)*t + e1)*t + e0) / x
  RETURN
END IF

!                        X >= 50

t = (1._dp/x) ** 2
z = (((((((((((s11*t + s10)*t + s9)*t + s8)*t + s7)*t + s6)*t + s5)*t +  &
    s4)*t + s3)*t + s2)*t + s1)*t - 0.5_dp) * t + 1._dp
fn_val = rpinv * (z/x)
RETURN
END FUNCTION derfc0



FUNCTION drexp(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!            EVALUATION OF THE FUNCTION EXP(X) - 1
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN)  :: x
REAL (dp)              :: fn_val

! Local variables
REAL (dp)  :: e, w, z
REAL (dp)  :: a0 = .248015873015873015873016D-04,   &
    a1 = -.344452080605731005808147D-05, a2 = .206664230430046597475413D-06,  &
    a3 = -.447300111094328162971036D-08, a4 = .114734027080634968083920D-11,  &
    b1 = -.249994190011341852652396_dp, b2 = .249987228833107957725728D-01,  &
    b3 = -.119037506846942249362528D-02, b4 = .228908693387350391768682D-04
REAL (dp) :: c1 = .1666666666666666666666666666666667_dp,   &
             c2 = .4166666666666666666666666666666667D-01,   &
             c3 = .8333333333333333333333333333333333D-02,   &
             c4 = .1388888888888888888888888888888889D-02,   &
             c5 = .1984126984126984126984126984126984D-03
!---------------------------
IF (ABS(x) <= 0.15_dp) THEN

!     Z IS A MINIMAX APPROXIMATION OF THE SERIES

!             C6 + C7*X + C8*X**2 + ....

!     THIS APPROXIMATION IS ACCURATE TO WITHIN
!     1 UNIT OF THE 23-RD SIGNIFICANT DIGIT.
!     THE RESULTING VALUE FOR W IS ACCURATE TO
!     WITHIN 1 UNIT OF THE 33-RD SIGNIFICANT DIGIT.

  z = ((((a4*x + a3)*x + a2)*x + a1)*x + a0) /  &
      ((((b4*x + b3)*x + b2)*x + b1)*x + 1._dp)
  w = ((((((z*x + c5)*x + c4)*x + c3)*x + c2)*x + c1)*x + 0.5_dp)*x + 1._dp
  fn_val = x * w
  RETURN
END IF

IF (x >= 0._dp) THEN
  e = EXP(x)
  fn_val = e * (0.5_dp + (0.5_dp - 1._dp/e))
  RETURN
END IF
IF (x >= -77._dp) THEN
  fn_val = (EXP(x) - 0.5_dp) - 0.5_dp
  RETURN
END IF
fn_val = -1._dp
RETURN
END FUNCTION drexp



FUNCTION dlnrel(a) RESULT(fn_val)
!-----------------------------------------------------------------------
!            EVALUATION OF THE FUNCTION LN(1 + A)
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN)  :: a
REAL (dp)              :: fn_val

! Local variables
REAL (dp) :: t, t2, w, z
REAL (dp) :: p0 = .7692307692307692307680D-01,   &
       p1 = -.1505958055914600184836_dp, p2 =  .9302355725278521726,   &
       p3 = -.1787900022182327735804D-01, q1 = -.2824412139355646910683D+01,  &
       q2 =  .2892424216041495392509D+01, q3 = -.1263560605948009364422D+01,  &
       q4 =  .1966769435894561313526_dp
REAL (dp) :: c1 = .3333333333333333333333333333333_dp,   &
             c2 = .2000000000000000000000000000000_dp,   &
             c3 = .1428571428571428571428571428571_dp,   &
             c4 = .1111111111111111111111111111111_dp,   &
             c5 = .9090909090909090909090909090909D-01
!-------------------------
IF (ABS(a) >= 0.375_dp) THEN
  t = 1._dp + a
  IF (a < 0._dp) t = 0.5_dp + (0.5_dp + a)
  fn_val = LOG(t)
  RETURN
END IF

!     W IS A MINIMAX APPROXIMATION OF THE SERIES

!            C6 + C7*T**2 + C8*T**4 + ...

!     THIS APPROXIMATION IS ACCURATE TO WITHIN 1.6 UNITS OF THE 21-ST
!     SIGNIFICANT DIGIT.
!     THE RESULTING VALUE FOR 1._dp + T2*Z IS ACCURATE TO WITHIN 1 UNIT OF
!     THE 30-TH SIGNIFICANT DIGIT.

t = a / (a + 2._dp)
t2 = t * t
w = (((p3*t2 + p2)*t2 + p1)*t2 + p0) /  &
    ((((q4*t2 + q3)*t2 + q2)*t2 + q1)*t2 + 1._dp)
z = ((((w*t2 + c5)*t2 + c4)*t2 + c3)*t2 + c2)*t2 + c1
fn_val = 2._dp * t * (1._dp + t2*z)
RETURN
END FUNCTION dlnrel


FUNCTION dxparg (l) RESULT(fn_val)
!--------------------------------------------------------------------
!     IF L = 0 THEN  DXPARG(L) = THE LARGEST POSITIVE W FOR WHICH
!     DEXP(W) CAN BE COMPUTED.
!
!     IF L IS NONZERO THEN  DXPARG(L) = THE LARGEST NEGATIVE W FOR
!     WHICH THE COMPUTED VALUE OF DEXP(W) IS NONZERO.
!
!     NOTE... ONLY AN APPROXIMATE VALUE FOR DXPARG(L) IS NEEDED.
!--------------------------------------------------------------------
IMPLICIT NONE
INTEGER, INTENT(IN) :: l
REAL (dp)           :: fn_val

! Local variable
REAL (dp)    :: one = 1._dp

IF (l == 0) THEN
  fn_val = LOG( HUGE(one) )
ELSE
  fn_val = LOG( TINY(one) )
END IF
RETURN
END FUNCTION dxparg


FUNCTION drlog(x) RESULT(fn_val)
!--------------------------------------------------------------------
!          EVALUATION OF THE FUNCTION X - 1 - LN(X)
!--------------------------------------------------------------------

REAL (dp), INTENT(IN)  :: x
REAL (dp)              :: fn_val

! Local variables
REAL (dp), PARAMETER :: a = .566749439387323789126387112411845D-01,   &
                        b = .456512608815524058941143273395059D-01
REAL (dp), PARAMETER :: p0 = .7692307692307692307680D-01,   &
       p1 = -.1505958055914600184836_dp, p2 = .9302355725278521726,   &
       p3 = -.1787900022182327735804D-01, q1 = -.2824412139355646910683D+01, &
       q2 = .2892424216041495392509D+01, q3 = -.1263560605948009364422D+01,  &
       q4 = .1966769435894561313526_dp
REAL (dp), PARAMETER :: c1 = .333333333333333333333333333333333_dp,   &
                        c2 = .200000000000000000000000000000000_dp,   &
                        c3 = .142857142857142857142857142857143_dp,   &
                        c4 = .111111111111111111111111111111111_dp,   &
                        c5 = .909090909090909090909090909090909D-01
REAL (dp)  :: r, t, u, up2, w, w1, z
!-------------------------
!     A = DRLOG (0.7)
!     B = DRLOG (4/3)
!-------------------------
IF (x >= 0.61_dp .AND. x <= 1.57_dp) THEN
  IF (x >= 0.82_dp) THEN
    IF (x > 1.18_dp) GO TO 10

!                 ARGUMENT REDUCTION

    u = (x-0.5_dp) - 0.5_dp
    up2 = u + 2._dp
    w1 = 0._dp
    GO TO 20
  END IF

  u = (x-0.7_dp) / 0.7_dp
  up2 = u + 2._dp
  w1 = a - u * 0.3_dp
  GO TO 20

  10 t = 0.75_dp * (x-1._dp)
  u = t - 0.25_dp
  up2 = t + 1.75_dp
  w1 = b + u / 3._dp

!                  SERIES EXPANSION

  20 r = u / up2
  t = r * r

!     Z IS A MINIMAX APPROXIMATION OF THE SERIES

!            C6 + C7*R**2 + C8*R**4 + ...

!     FOR THE INTERVAL (0.0, 0.375). THE APPROXIMATION IS ACCURATE
!     TO WITHIN 1.6 UNITS OF THE 21-ST SIGNIFICANT DIGIT.

  z = (((p3*t + p2)*t + p1)*t + p0) / ((((q4*t + q3)*t + q2)*t + q1)*t + 1._dp)

  w = ((((z*t + c5)*t + c4)*t + c3)*t + c2) * t + c1
  fn_val = r * (u-2._dp*t*w) + w1
  RETURN
END IF

r = (x-0.5_dp) - 0.5_dp
fn_val = r - LOG(x)
RETURN
END FUNCTION drlog



FUNCTION dpdel(x) RESULT(fn_val)
!--------------------------------------------------------------------

!  COMPUTATION OF THE FUNCTION DEL(X) FOR  X >= 10  WHERE
!  LN(GAMMA(X)) = (X - 0.5)*LN(X) - X + 0.5*LN(2*PI) + DEL(X)

!                      --------

!  THE SERIES FOR DPDEL ON THE INTERVAL 0.0 TO 1.0 DERIVED BY
!  A.H. MORRIS FROM THE CHEBYSHEV SERIES IN THE SLATEC LIBRARY
!  OBTAINED BY WAYNE FULLERTON (LOS ALAMOS).

!--------------------------------------------------------------------
REAL (dp), INTENT(IN)  :: x
REAL (dp)              :: fn_val

! Local variables
REAL (dp), PARAMETER :: a(15) = (/ .833333333333333333333333333333D-01,  &
        -.277777777777777777777777752282D-04,  &
         .793650793650793650791732130419D-07,  &
        -.595238095238095232389839236182D-09,  &
         .841750841750832853294451671990D-11,  &
        -.191752691751854612334149171243D-12,  &
         .641025640510325475730918472625D-14,  &
        -.295506514125338232839867823991D-15,  &
         .179643716359402238723287696452D-16,  &
        -.139228964661627791231203060395D-17,  &
         .133802855014020915603275339093D-18,  &
        -.154246009867966094273710216533D-19,  &
         .197701992980957427278370133333D-20,  &
        -.234065664793997056856992426667D-21,  &
         .171348014966398575409015466667D-22 /)
REAL (dp) :: t, w
INTEGER   :: i, k
!-----------------------------------------------------------------------
t = (10._dp/x) ** 2
w = a(15)
DO i = 1, 14
  k = 15 - i
  w = t * w + a(k)
END DO
fn_val = w / x
RETURN
END FUNCTION dpdel



FUNCTION dsin1(x) RESULT(fn_val)
!--------------------------------------------------------------------

!             REAL (dp) EVALUATION OF SIN(PI*X)

!                          --------------

!  THE EXPANSION FOR SIN(PI*A) (ABS(A) <= PI/4) USING A1,...,A13
!  IS ACCURATE TO WITHIN 2 UNITS OF THE 40-TH SIGNIFICANT DIGIT, AND
!  THE EXPANSION FOR COS(PI*A) (ABS(A) <= PI/4) USING B1,...,B13
!  IS ACCURATE TO WITHIN 4 UNITS OF THE 40-TH SIGNIFICANT DIGIT.

!--------------------------------------------------------------------
REAL (dp), INTENT(IN)  :: x
REAL (dp)              :: fn_val

! Local variables
REAL (dp)            :: a, t, w
REAL (dp), PARAMETER :: pi = 3.141592653589793238462643383279502884197_dp
REAL (dp), PARAMETER :: a1 = -.1028083791780141522795259479153765743002_dp,   &
      a2  = .3170868848763100170457042079710451905600D-02,   &
      a3  = -.4657026956105571623449026167864697920000D-04,  &
      a4  = .3989844942879455643410226655783424000000D-06,   &
      a5  = -.2237397227721999776371894030796800000000D-08,  &
      a6  = .8847045483056962709715066675200000000000D-11,   &
      a7  = -.2598715447506450292885585920000000000000D-13,  &
      a8  = .5893449774331011070033920000000000000000D-16 ,  &
      a9  = -.1062975472045522550784000000000000000000D-18,   &
      a10 = .1561182648301780992000000000000000000000D-21,    &
      a11 = -.1903193516670976000000000000000000000000D-24,   &
      a12 = .1956617650176000000000000000000000000000D-27,    &
      a13 = -.1711276032000000000000000000000000000000D-30
REAL (dp), PARAMETER :: b1 = -.3084251375340424568385778437461297229882_dp, &
      b2  = .1585434424381550085228521039855226435920D-01,   &
      b3  = -.3259918869273900136414318317506279360000D-03,  &
      b4  = .3590860448591510079069203991239232000000D-05,   &
      b5  = -.2461136950494199754009084061808640000000D-07,  &
      b6  = .1150115912797405152263195572224000000000D-09,   &
      b7  = -.3898073171259675439899172864000000000000D-12,  &
      b8  = .1001886461636271969091584000000000000000D-14,   &
      b9  = -.2019653396886572027084800000000000000000D-17,  &
      b10 = .3278483561466560512000000000000000000000D-20,   &
      b11 = -.4377345082051788800000000000000000000000D-23,  &
      b12 = .4891532381388800000000000000000000000000D-26,   &
      b13 = -.4617089843200000000000000000000000000000D-29
INTEGER  :: max, n
!------------------------

!     ****** MAX IS A MACHINE DEPENDENT CONSTANT. MAX IS THE
!            LARGEST POSITIVE INTEGER THAT MAY BE USED.

!                       MAX = IPMPAR(3)
max = HUGE(3)

!------------------------
a = ABS(x)
t = MAX
IF (a >= t) THEN
  fn_val = 0._dp
  RETURN
END IF

n = a
t = n
a = a - t
IF (a <= 0.75_dp) THEN
  IF (a < 0.25_dp) GO TO 10

!                    0.25 <= A <= 0.75

  a = 0.25_dp + (0.25_dp-a)
  t = 16._dp * a * a
  fn_val = (((((((((((((b13*t + b12)*t + b11)*t + b10)*t + b9)*t + b8)*t  &
           + b7)*t + b6)*t + b5)*t + b4)*t + b3)*t + b2)*t + b1)*t +  &
           0.5_dp) + 0.5_dp
  GO TO 20
END IF

!                 A < 0.25  OR  A > 0.75

a = 0.25_dp + (0.75_dp-a)
10 t = 16._dp * a * a
w = (((((((((((((a13*t + a12)*t + a11)*t + a10)*t + a9)*t + a8)*t + a7)*t  &
    + a6)*t + a5)*t + a4)*t + a3)*t + a2)*t + a1)*t + 0.5_dp) + 0.5_dp
fn_val = pi * a * w

!                        TERMINATION

20 IF (x < 0.0) fn_val = -fn_val
IF (MOD(n,2) /= 0) fn_val = -fn_val
RETURN
END FUNCTION dsin1



FUNCTION dgamma(a) RESULT(fn_val)
!--------------------------------------------------------------------

!             EVALUATION OF THE GAMMA FUNCTION FOR
!                  REAL (dp) ARGUMENTS

!                        -----------

!  DGAMMA(A) IS ASSIGNED THE VALUE 0 WHEN THE GAMMA FUNCTION CANNOT
!  BE COMPUTED.

!--------------------------------------------------------------------
!  WRITTEN BY ALFRED H. MORRIS, JR.
!       NAVAL SURFACE WEAPONS CENTER
!       DAHLGREN, VIRGINIA
!--------------------------------------------------------------------
REAL (dp), INTENT(IN)  :: a
REAL (dp)              :: fn_val

! Local variables
REAL (dp), PARAMETER :: d = 0.41893853320467274178032973640562_dp,  &
                        pi = 3.14159265358979323846264338327950_dp
REAL (dp) :: s, t, x, w
INTEGER   :: j, n
!-----------------------------------------------------------------------
!     D = 0.5*(LN(2*PI) - 1)
!-----------------------------------------------------------------------
fn_val = 0._dp
x = a
IF (ABS(a) <= 20._dp) THEN
!-----------------------------------------------------------------------
!             EVALUATION OF DGAMMA(A) FOR ABS(A) <= 20
!-----------------------------------------------------------------------
  t = 1._dp
  n = x
  n = n - 1

!     LET T BE THE PRODUCT OF A-J WHEN A >= 2

  IF (n < 0) THEN
    GO TO 40
  ELSE IF (n == 0) THEN
    GO TO 30
  END IF

  DO j = 1, n
    x = x - 1._dp
    t = x * t
  END DO
  30 x = x - 1._dp
  GO TO 60

!     LET T BE THE PRODUCT OF A+J WHEN A < 1

  40 t = a
  IF (a <= 0._dp) THEN
    n = -n - 1
    IF (n /= 0) THEN
      DO j = 1, n
        x = x + 1._dp
        t = x * t
      END DO
    END IF
    x = (x+0.5_dp) + 0.5_dp
    t = x * t
    IF (t == 0._dp) RETURN
  END IF


!     THE FOLLOWING CODE CHECKS IF 1/T CAN OVERFLOW. THIS
!     CODE MAY BE OMITTED IF DESIRED.

  IF (ABS(t) < 1.d-33) THEN
    IF (ABS(t)*HUGE(1.0_dp) <= 1.000000001_dp) RETURN
    fn_val = 1._dp / t
    RETURN
  END IF

!     COMPUTE DGAMMA(1 + X) FOR 0 <= X < 1

  60 fn_val = 1._dp / (1._dp + dgam1(x))

!     TERMINATION

  IF (a >= 1._dp) THEN
    fn_val = fn_val * t
    RETURN
  END IF
  fn_val = fn_val / t
  RETURN
END IF
!-----------------------------------------------------------------------
!           EVALUATION OF DGAMMA(A) FOR ABS(A) > 20
!-----------------------------------------------------------------------
IF (ABS(a) >= 1.d3) RETURN
IF (a <= 0._dp) THEN
  s = dsin1(a) / pi
  IF (s == 0._dp) RETURN
  x = -a
END IF

!     COMPUTE THE MODIFIED ASYMPTOTIC SUM

w = dpdel(x)

!     FINAL ASSEMBLY

w = (d+w) + (x-0.5_dp) * (LOG(x)-1._dp)
IF (w > dxparg(0)) RETURN
fn_val = EXP(w)
IF (a < 0._dp) fn_val = (1._dp/(fn_val*s)) / x

RETURN
END FUNCTION dgamma



FUNCTION dgam1(x) RESULT(fn_val)
!--------------------------------------------------------------------
!  EVALUATION OF 1/GAMMA(1 + X) - 1  FOR -0.5 <= X <= 1.5
!--------------------------------------------------------------------

!  THE FOLLOWING ARE THE FIRST 49 COEFFICIENTS OF THE MACLAURIN
!  EXPANSION FOR 1/GAMMA(1 + X) - 1. THE COEFFICIENTS ARE
!  CORRECT TO 40 DIGITS.  THE COEFFICIENTS WERE OBTAINED BY
!  ALFRED H. MORRIS JR. (NAVAL SURFACE WARFARE CENTER) AND ARE
!  GIVEN HERE FOR REFERENCE.  ONLY THE FIRST 14 COEFFICIENTS ARE
!  USED IN THIS CODE.

!                        -----------

!  DATA A(1)  / .5772156649015328606065120900824024310422_dp/,
! *     A(2)  /-.6558780715202538810770195151453904812798_dp/,
! *     A(3)  /-.4200263503409523552900393487542981871139D-01/,
! *     A(4)  / .1665386113822914895017007951021052357178_dp/,
! *     A(5)  /-.4219773455554433674820830128918739130165D-01/,
! *     A(6)  /-.9621971527876973562114921672348198975363D-02/,
! *     A(7)  / .7218943246663099542395010340446572709905D-02/,
! *     A(8)  /-.1165167591859065112113971084018388666809D-02/,
! *     A(9)  /-.2152416741149509728157299630536478064782D-03/,
! *     A(10) / .1280502823881161861531986263281643233949D-03/
!  DATA A(11) /-.2013485478078823865568939142102181838229D-04/,
! *     A(12) /-.1250493482142670657345359473833092242323D-05/,
! *     A(13) / .1133027231981695882374129620330744943324D-05/,
! *     A(14) /-.2056338416977607103450154130020572836513D-06/,
! *     A(15) / .6116095104481415817862498682855342867276D-08/,
! *     A(16) / .5002007644469222930055665048059991303045D-08/,
! *     A(17) /-.1181274570487020144588126565436505577739D-08/,
! *     A(18) / .1043426711691100510491540332312250191401D-09/,
! *     A(19) / .7782263439905071254049937311360777226068D-11/,
! *     A(20) /-.3696805618642205708187815878085766236571D-11/
!  DATA A(21) / .5100370287454475979015481322863231802727D-12/,
! *     A(22) /-.2058326053566506783222429544855237419746D-13/,
! *     A(23) /-.5348122539423017982370017318727939948990D-14/,
! *     A(24) / .1226778628238260790158893846622422428165D-14/,
! *     A(25) /-.1181259301697458769513764586842297831212D-15/,
! *     A(26) / .1186692254751600332579777242928674071088D-17/,
! *     A(27) / .1412380655318031781555803947566709037086D-17/,
! *     A(28) /-.2298745684435370206592478580633699260285D-18/,
! *     A(29) / .1714406321927337433383963370267257066813D-19/,
! *     A(30) / .1337351730493693114864781395122268022875D-21/
!  DATA A(31) /-.2054233551766672789325025351355733796682D-21/,
! *     A(32) / .2736030048607999844831509904330982014865D-22/,
! *     A(33) /-.1732356445910516639057428451564779799070D-23/,
! *     A(34) /-.2360619024499287287343450735427531007926D-25/,
! *     A(35) / .1864982941717294430718413161878666898946D-25/,
! *     A(36) /-.2218095624207197204399716913626860379732D-26/,
! *     A(37) / .1297781974947993668824414486330594165619D-27/,
! *     A(38) / .1180697474966528406222745415509971518560D-29/,
! *     A(39) /-.1124584349277088090293654674261439512119D-29/,
! *     A(40) / .1277085175140866203990206677751124647749D-30/
!  DATA A(41) /-.7391451169615140823461289330108552823711D-32/,
! *     A(42) / .1134750257554215760954165259469306393009D-34/,
! *     A(43) / .4639134641058722029944804907952228463058D-34/,
! *     A(44) /-.5347336818439198875077418196709893320905D-35/,
! *     A(45) / .3207995923613352622861237279082794391090D-36/,
! *     A(46) /-.4445829736550756882101590352124643637401D-38/,
! *     A(47) /-.1311174518881988712901058494389922190237D-38/,
! *     A(48) / .1647033352543813886818259327906394145400D-39/,
! *     A(49) /-.1056233178503581218600561071538285049997D-40/

!                        -----------

!  C = A(1) - 1 IS ALSO FREQUENTLY NEEDED. C HAS THE VALUE ...

!  DATA C /-.4227843350984671393934879099175975689578_dp/

!--------------------------------------------------------------------
REAL (dp), INTENT(IN) :: x
REAL (dp)             :: fn_val

! Local variables
REAL (dp) :: d, t, w, z
REAL (dp), PARAMETER :: a0 = .611609510448141581788D-08, a1  &
        = .624730830116465516210D-08, b1 = .203610414066806987300_dp, b2  &
        = .266205348428949217746D-01, b3 = .493944979382446875238D-03, b4  &
        = -.851419432440314906588D-05, b5 = -.643045481779353022248D-05, b6  &
        = .992641840672773722196D-06, b7 = -.607761895722825260739D-07, b8  &
        = .195755836614639731882D-09
REAL (dp), PARAMETER :: p0 = .6116095104481415817861D-08, p1  &
        = .6871674113067198736152D-08, p2 = .6820161668496170657, p3  &
        = .4686843322948848031080D-10, p4 = .1572833027710446286995D-11, p5  &
        = -.1249441572276366213222D-12, p6 = .4343529937408594255178D-14, q1  &
        = .3056961078365221025009_dp, q2 = .5464213086042296536016D-01, q3  &
        = .4956830093825887312, q4 = .2692369466186361192876D-03
REAL (dp), PARAMETER :: c = -.422784335098467139393487909917598_dp, c0  &
        = .577215664901532860606512090082402_dp, c1  &
        = -.655878071520253881077019515145390_dp, c2  &
        = -.420026350340952355290039348754298D-01, c3  &
        = .166538611382291489501700795102105_dp, c4  &
        = -.421977345555443367482083012891874D-01, c5  &
        = -.962197152787697356211492167234820D-02, c6  &
        = .721894324666309954239501034044657D-02, c7  &
        = -.116516759185906511211397108401839D-02, c8  &
        = -.215241674114950972815729963053648D-03, c9  &
        = .128050282388116186153198626328164D-03, c10  &
        = -.201348547807882386556893914210218D-04, c11  &
        = -.125049348214267065734535947383309D-05, c12  &
        = .113302723198169588237412962033074D-05, c13  &
        = -.205633841697760710345015413002057D-06
!----------------------------
t = x
d = x - 0.5_dp
IF (d > 0._dp) t = d - 0.5_dp
IF (t < 0.0_dp) THEN
  GO TO 30
ELSE IF (t > 0.0_dp) THEN
  GO TO 20
END IF

fn_val = 0._dp
RETURN
!------------

!             CASE WHEN 0 < T <= 0.5

!           W IS A MINIMAX APPROXIMATION FOR
!           THE SERIES A(15) + A(16)*T + ...

!------------
20 w = ((((((p6*t + p5)*t + p4)*t + p3)*t + p2)*t + p1)*t + p0) /   &
       ((((q4*t+q3)*t + q2)*t + q1)*t + 1._dp)
z = (((((((((((((w*t + c13)*t + c12)*t + c11)*t + c10)*t + c9)*t + c8)*t + c7)*t  &
    + c6)*t + c5)*t + c4)*t + c3)*t + c2)*t + c1) * t + c0

IF (d <= 0._dp) THEN
  fn_val = x * z
  RETURN
END IF
fn_val = (t/x) * ((z-0.5_dp)-0.5_dp)
RETURN
!------------

!             CASE WHEN -0.5 <= T < 0

!           W IS A MINIMAX APPROXIMATION FOR
!           THE SERIES A(15) + A(16)*T + ...

!------------
30 w = (a1*t + a0) / ((((((((b8*t + b7)*t + b6)*t + b5)*t + b4)*t + b3)*t + b2)*t + b1)*t+1._dp)
z = (((((((((((((w*t + c13)*t + c12)*t + c11)*t + c10)*t + c9)*t + c8)*t + c7)*t  &
    + c6)*t + c5)*t + c4)*t + c3)*t + c2)*t + c1) * t + c

IF (d <= 0._dp) THEN
  fn_val = x * ((z+0.5_dp)+0.5_dp)
  RETURN
END IF
fn_val = t * z / x
RETURN
END FUNCTION dgam1



FUNCTION dgamln(a) RESULT(fn_val)
!--------------------------------------------------------------------

!        EVALUATION OF LN(GAMMA(A)) FOR POSITIVE A

!--------------------------------------------------------------------
!  WRITTEN BY ALFRED H. MORRIS
!       NAVAL SURFACE WEAPONS CENTER
!       DAHLGREN, VIRGINIA
!--------------------------------------------------------------------
!  D = 0.5*(LN(2*PI) - 1)
!--------------------------------------------------------------------
REAL (dp), INTENT(IN)  :: a
REAL (dp)              :: fn_val

! Local variables
REAL (dp), PARAMETER  :: d = 0.41893853320467274178032973640562_dp
REAL (dp) :: w, x
INTEGER   :: i, n
!--------------------------
IF (a < 0.5_dp) THEN
  fn_val = dgmln1(a) - LOG(a)
  RETURN
END IF
IF (a <= 2.5_dp) THEN
  x = a - 1._dp
  IF (a < 1._dp) x = (a-0.5_dp) - 0.5_dp
  fn_val = dgmln1(x)
  RETURN
END IF

IF (a < 10._dp) THEN
  n = a - 1.5_dp
  x = a
  w = 1._dp
  DO i = 1, n
    x = x - 1._dp
    w = x * w
  END DO
  fn_val = dgmln1(x-1._dp) + LOG(w)
  RETURN
END IF

w = dpdel(a)
fn_val = (d+w) + (a-0.5_dp) * (LOG(a)-1._dp)
RETURN
END FUNCTION dgamln



FUNCTION dgmln1(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!     EVALUATION OF LN(GAMMA(1 + X)) FOR -0.5 <= X <= 1.5
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN)  :: x
REAL (dp)              :: fn_val

! Local variables
REAL (dp)  :: w
!-----------------------
w = dgam1(x)
fn_val = -dlnrel(w)
RETURN
END FUNCTION dgmln1



FUNCTION drcomp(a, x) RESULT(fn_val)
!-----------------------------------------------------------------------
!              EVALUATION OF EXP(-X)*X**A/GAMMA(A)
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN)  :: a, x
REAL (dp)              :: fn_val

! Local variables
REAL (dp), PARAMETER  :: c = .398942280401432677939946059934_dp
REAL (dp)             :: t, w
!--------------------------
!     C = 1/SQRT(2*PI)
!--------------------------
fn_val = 0._dp
IF (x == 0._dp) RETURN
IF (a <= 20._dp) THEN

  t = a * LOG(x) - x
  IF (t < dxparg(1)) RETURN
  IF (a < 1._dp) THEN
    fn_val = (a*EXP(t)) * (1._dp + dgam1(a))
    RETURN
  END IF
  fn_val = EXP(t) / dgamma(a)
  RETURN
END IF

t = x / a
IF (t == 0._dp) RETURN
w = -(dpdel(a) + a*drlog(t))
IF (w >= dxparg(1)) fn_val = c * SQRT(a) * EXP(w)
RETURN
END FUNCTION drcomp



SUBROUTINE dpni(p, q, d, w, ierr)
!--------------------------------------------------------------------

!      EVALUATION OF THE INVERSE NORMAL DISTRIBUTION FUNCTION

!                        ------------

!  LET F(T) = 1/(SQRT(2*PI)*EXP(-T*T/2)). THEN THE FUNCTION

!     PROB(X) = INTEGRAL FROM MINUS INFINITY TO X OF F(T)

!  IS THE NORMAL DISTRIBUTION FUNCTION OF ZERO MEAN AND UNIT
!  VARIANCE. IT IS ASSUMED THAT P > 0, Q > 0, P + Q = 1,
!  AND D = P - 0.5. THE VALUE W IS COMPUTED WHERE PROB(W) = P.

!  IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.

!    IERR = 0  NO INPUT ERRORS WERE DETECTED. W WAS COMPUTED.
!    IERR = 1  EITHER P OR Q IS INCORRECT.
!    IERR = 2  D IS INCORRECT.

!--------------------------------------------------------------------
REAL (dp), INTENT(IN)   :: p, q, d
REAL (dp), INTENT(OUT)  :: w
INTEGER, INTENT(OUT)    :: ierr

! Local variables
REAL (dp), PARAMETER  :: rt2 = 1.4142135623730950488016887242097_dp
REAL (dp)             :: eps, t, u, v
!------------------------
!     RT2 = SQRT(2)
!------------------------
t = MIN(p,q)
IF (t > 0._dp) THEN
  eps = EPSILON(1.0_dp)
  w = 0.5_dp + (0.5_dp - (p+q))
  IF (ABS(w) <= 2._dp*eps) THEN

    u = ABS(d+d)
    v = t + t
    w = derfi(u,v)
    IF (w < 0._dp) GO TO 10

    ierr = 0
    w = rt2 * w
    IF (d < 0._dp) w = -w
    RETURN
  END IF
END IF

!                         ERROR RETURN

ierr = 1
RETURN
10 ierr = 2
RETURN
END SUBROUTINE dpni



FUNCTION derfi(p, q) RESULT(fn_val)
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN)  :: p, q
REAL (dp)              :: fn_val

!               REAL (dp) COMPUTATION OF
!                 THE INVERSE ERROR FUNCTION

!                      ----------------

!  FOR 0 <= P <= 1,  W = DERFI(P,Q) WHERE ERF(W) = P. IT
!  IS ASSUMED THAT Q = 1 - P. IF P < 0, Q <= 0, OR P + Q
!  IS NOT 1, THEN DERFI(P,Q) IS SET TO A NEGATIVE VALUE.

!--------------------------------------------------------------------
!  REFERENCE. MATHEMATICS OF COMPUTATION,OCT.1976,PP.827-830.
!               J.M.BLAIR,C.A.EDWARDS,J.H.JOHNSON
!--------------------------------------------------------------------
REAL (dp) :: c = .5625_dp, c1 = .87890625_dp, c2  &
        = -.2302585092994045684017991454684364D+03, r  &
        = .8862269254527580136490837416705726D+00, eps, f, lnq, s, t, x
REAL (dp) :: a(7) = (/ .841467547194693616D-01,  &
        .160499904248262200D+01, .809451641478547505D+01,  &
        .164273396973002581D+02, .154297507839223692D+02,  &
        .669584134660994039D+01, .108455979679682472D+01 /), a1(7)  &
        = (/ .552755110179178015D+2, .657347545992519152D+3,  &
        .124276851197202733D+4, .818859792456464820D+3,  &
        .234425632359410093D+3, .299942187305427917D+2,  &
        .140496035731853946D+1 /), a2(7) = (/ .500926197430588206D+1,  &
        .111349802614499199D+3, .353872732756132161D+3,  &
        .356000407341490731D+3, .143264457509959760D+3,  &
        .240823237485307567D+2, .140496035273226366D+1 /), a3(11)  &
        = (/ .237121026548776092D4, .732899958728969905D6,  &
        .182063754893444775D7, .269191299062422172D7, .304817224671614253D7  &
       , .130643103351072345D7, .296799076241952125D6,  &
        .457006532030955554D5, .373449801680687213D4, .118062255483596543D3  &
       , .100000329157954960D1 /), a4(9) = (/ .154269429680540807D12,  &
        .430207405012067454D12, .182623446525965017D12,  &
        .248740194409838713D11, .133506080294978121D10,  &
        .302446226073105850D08, .285909602878724425D06,  &
        .101789226017835707D04, .100000004821118676D01 /),  &
        b(7) = (/ .352281538790042405D-02, .293409069065309557D+00,  &
        .326709873508963100D+01, .123611641257633210D+02,  &
        .207984023857547070D+02, .170791197367677668D+02,  &
        .669253523595376683D+01 /), b1(6) = (/ .179209835890172156D+3,  &
        .991315839349539886D+3, .138271033653003487D+4, .764020340925985926D+3,  &
        .194354053300991923D+3, .228139510050586581D+2 /),  &
        b2(6) = (/ .209004294324106981D+2, .198607335199741185D+3,  &
        .439311287748524270D+3, .355415991280861051D+3, .123303672628828521D+3,  &
        .186060775181898848D+2 /), b3(10) = (/ .851911109952055378D6,  &
        .194746720192729966D7, .373640079258593694D7, .397271370110424145D7,  &
        .339457682064283712D7 , .136888294898155938D7, .303357770911491406D6,  &
        .459721480357533823D5, .373762573565814355D4, .118064334590001264D3 /),  &
        b4(9) = (/ .220533001293836387D12, .347822938010402687D12,  &
        .468373326975152250D12, .185251723580351631D12,  &
        .249464490520921771D11, .133587491840784926D10,  &
        .302480682561295591D08, .285913799407861384D06,  &
        .101789250893050230D04 /)
!-----------------------------------------------------------------------
!     C2 = LN(1.E-100)
!     R  = SQRT(PI)/2
!-----------------------------------------------------------------------
IF (p >= 0._dp .AND. q > 0._dp) THEN
  eps = EPSILON(1.0_dp)
  t = 0.5_dp + (0.5_dp-(p+q))
  IF (ABS(t) > 3._dp*eps) GO TO 10

!                      0 <= P <= 0.75

  IF (p <= 0.75_dp) THEN
    x = c - p * p
    s = (((((a(1)*x+a(2))*x+a(3))*x+a(4))*x+a(5))*x+a(6)) * x +a(7)
    t = ((((((b(1)*x+b(2))*x+b(3))*x+b(4))*x+b(5))*x+b(6))*x+b(7)) * x + 1._dp
    fn_val = p * (s/t)
    IF (eps > 1.d-19) RETURN

    x = fn_val
    f = derf(x) - p
    fn_val = x - r * EXP(x*x) * f
    RETURN
  END IF

!                    0.75 < P <= 0.9375

  IF (p <= 0.9375_dp) THEN
    x = c1 - p * p
    IF (x <= 0.1_dp) THEN
      s = ((((((a1(1)*x+a1(2))*x+a1(3))*x+a1(4))*x+a1(5))*x+a1(6))*x+a1(7))
      t = ((((((b1(1)*x+b1(2))*x+b1(3))*x+b1(4))*x+b1(5))*x+b1(6))*x+1._dp)
    ELSE

      s = ((((((a2(1)*x+a2(2))*x+a2(3))*x+a2(4))*x+a2(5))*x+a2(6))*x+a2(7))
      t = ((((((b2(1)*x+b2(2))*x+b2(3))*x+b2(4))*x+b2(5))*x+b2(6))*x+1._dp)
    END IF

    fn_val = p * (s/t)
    IF (eps > 1.d-19) RETURN

    x = fn_val
    t = derfc1(1,x) - EXP(x*x) * q
    fn_val = x + r * t
    RETURN
  END IF

!                  1.E-100 <= Q < 0.0625

  lnq = LOG(q)
  x = 1._dp / SQRT(-lnq)
  IF (lnq >= c2) THEN
    s = (((((((((a3(1)*x+a3(2))*x+a3(3))*x+a3(4))*x+a3(5))*x+  &
    a3(6))*x+a3(7))*x+a3(8))*x+a3(9))*x+a3(10)) * x + a3(11)
    t = (((((((((b3(1)*x+b3(2))*x+b3(3))*x+b3(4))*x+b3(5))*x+  &
    b3(6))*x+b3(7))*x+b3(8))*x+b3(9))*x+b3(10)) * x + 1._dp
  ELSE

!                 1.E-10000 <= Q < 1.E-100

    s = (((((((a4(1)*x+a4(2))*x+a4(3))*x+a4(4))*x+a4(5))*x+a4(6))*  &
    x+a4(7))*x+a4(8)) * x + a4(9)
    t = ((((((((b4(1)*x+b4(2))*x+b4(3))*x+b4(4))*x+b4(5))*x+  &
    b4(6))*x+b4(7))*x+b4(8))*x+b4(9)) * x + 1._dp
  END IF

  fn_val = s / (x*t)
  IF (eps > 5.d-20) RETURN

  x = fn_val
  t = derfc1(1,x)
  f = (LOG(t)-lnq) - x * x
  fn_val = x + r * t * f
  RETURN
END IF

!                         ERROR RETURN

fn_val = -1._dp
RETURN
10 fn_val = -2._dp
RETURN
END FUNCTION derfi

END MODULE Incomplete_Gamma


!!****m* TOOLS/m_qtlmap_math
!!  NAME
!!    m_qtlmap_math -- Mathematical subroutines
!!  SYNOPSIS
!!    This modules purposes an interface with NAG specification. In a NAG environment, the nag subroutines is called
!!    otherwise an implementation of the subroutine using SLATEC is called
!!  DESCRIPTION
!!
!!  NOTES
!!
!!  BUGS
!!
!!  SEE ALSO
!!
!!  COPYRIGHT
!!***
!! You can use this space for remarks that should not be included
!! in the documentation.
!!/
module m_qtlmap_math
    use m_qtlmap_constant
    use m_qtlmap_log
    use m_qtlmap_base

    implicit none
    private

#ifndef HAVE_LIBNAG
    double precision    ,  parameter , private    :: SQRT_2 = 1.414213562373095d0 ! sqrt(2.d0)
#endif

    ! Give some information about library used by this module
    public :: MATH_QTLMAP_INFO
    ! computes the upper and lower tail probabilities and the probability density function of the beta
    ! distribution with parameters a and b.
    public :: MATH_QTLMAP_G01EEF

    ! calculates the mean, standard deviation, coefficients of skewness and kurtosis, and the maximum
    ! and minimum values for a set of ungrouped data. Weighting may be used.
    public :: MATH_QTLMAP_G01AAF

    !G01EAF returns a one or two-tail probability for the standard Normal distribution, via the routine name.
    public :: MATH_QTLMAP_G01EAF

    !G01FAF returns the deviate associated with the given probability of the standard Normal distribution, via
    ! the routine name.
    public :: MATH_QTLMAP_G01FAF

    ! returns the lower tail, upper tail or two-tail probability for the Student’s t-distribution with real
    ! degrees of freedom, via the routine name.
    public :: MATH_QTLMAP_G01EBF

    ! Performs a canonical variate (canonical discrimination) analysis.
    public :: MATH_QTLMAP_G03ACF

    public :: test_G03ACF

    ! performs a pseudo-random permutation of a vector of integers.
    public :: MATH_QTLMAP_G05EHF

    ! rearranges a vector of real numbers into ascending or descending order.
    public :: MATH_QTLMAP_M01CAF

    ! ranks a vector of real numbers in ascending or descending order.
    public :: MATH_QTLMAP_M01DAF

    ! Ranks a vector of integer numbers in ascending or descending order.
    public :: MATH_QTLMAP_M01DBF

    ! Ranks a vector of character data in ASCII or reverse ASCII order of a specified substring.
    public :: MATH_QTLMAP_M01DCF

    ! Inverts a permutation, and hence converts a rank vector to an index vector, or vice versa.
    public :: MATH_QTLMAP_M01ZAF

    ! Ranks the rows of a matrix of integer numbers in ascending or descending order.
    public :: MATH_QTLMAP_M01DFF

    ! returns the lower or upper tail probability for the 12 distribution with real degrees of freedom, via
    ! the routine name
    public :: MATH_QTLMAP_G01ECF

    !calculates the determinant of a real symmetric positive-definite matrix using a Cholesky factorization.
    public :: MATH_QTLMAP_F03ABF
    !calculates the inverse and the determinant of a real symmetric positive-definite.
    public :: MATH_QTLMAP_INVDETMATSYM
    !calculates the inverse and the determinant of a real matrix.
  !  public :: MATH_QTLMAP_INVDETMAT
    ! calculates the approximate inverse of a real symmetric positive-definite matrix, using a Cholesky
    ! factorization.
    public :: MATH_QTLMAP_F01ADF
    ! Computes the lower of the beta distribution with parameters a and b.
    public :: MATH_QTLMAP_LOWERTAIL_BETA

    !S15AEF returns the value of the error function erf x, via the routine name.
    public :: MATH_QTLMAP_S15AEF
    ! Returns the value of the complementary error function, erfc x, via the routine name.
    public :: MATH_QTLMAP_S15ADF
    !returns generalized inverse of matrix x of size n x n declared
    ! as m x m.
    public :: ginv1

    contains

!****************************************************************************************************************
!!****f* m_qtlmap_math/MATH_QTLMAP_INFO
!!  NAME
!!    MATH_QTLMAP_INFO
!!  DESCRIPTION
!!   Give some information about library used by this module
!!  NOTES
!!***
      subroutine MATH_QTLMAP_INFO

#ifdef HAVE_LIBNAG
           call log_mess("************* QTLMap use NAG library *****************************",INFO_DEF)
#else
           call log_mess("************* QTLMap use LAPACK/BLAS/Orderpack-2.0 library:*****************************",INFO_DEF)
#endif

      end subroutine MATH_QTLMAP_INFO
!****************************************************************************************************************
!!****f* m_qtlmap_math/MATH_QTLMAP_G01EEF
!!  NAME
!!    MATH_QTLMAP_G01EEF
!!  DESCRIPTION
!!   Computes the upper and lower tail probabilities and the probability density function of the beta
!!   distribution with parameters a and b.
!!  NOTES
!!***
      subroutine MATH_QTLMAP_G01EEF(X, A, B, TOL, P, Q, PDF, IFAIL)
           double precision,intent(in)     ::  X, A, B, TOL
           double precision,intent(out)    ::  P, Q, PDF
           integer , intent(inout)         ::  IFAIL

#ifndef HAVE_LIBNAG
           double precision                :: AGAMMA,BGAMMA
#endif

           call log_mess('MATH_QTLMAP_G01EEF (IN) ==== X:'//str(X),DEBUG_DEF)
           call log_mess('MATH_QTLMAP_G01EEF (IN) ==== A:'//str(A),DEBUG_DEF)
           call log_mess('MATH_QTLMAP_G01EEF (IN) ==== B:'//str(B),DEBUG_DEF)
           P=0.d0
           Q=0.d0
           PDF=0.d0

           IF ( X < 0 .or. X > 1 ) then
              IFAIL = 1
              CALL log_mess('MATH_QTLMAP_G01EEF:'//str(X)//' have to be between 0 and 1.',ERROR_DEF)
              RETURN
           END IF

           IF ( A < 0 )  then
             IFAIL = 2
             CALL log_mess('MATH_QTLMAP_G01EEF:A is less than 0',ERROR_DEF)
             RETURN
           END IF

           IF ( B < 0 ) then
             IFAIL = 2
             CALL log_mess('MATH_QTLMAP_G01EEF:B is less than 0',ERROR_DEF)
             RETURN
           END IF

#ifdef HAVE_LIBNAG
           call G01EEF(X, A, B, TOL, P, Q, PDF, IFAIL)
#else
           !The lower tail is also knwon as incomplete beta function
           P = DBETAI(X,A,B)

           IF (IFAIL /= 0 ) RETURN

           Q = 1 - P
           IF ( P == 0 .or. P == 1 ) THEN
             PDF = 0.d0
           ELSE
             AGAMMA = DGAMMA(A)
             BGAMMA = DGAMMA(B)
             IF (AGAMMA == 0 ) then
               IFAIL = 3
               CALL log_mess('MATH_QTLMAP_G01EEF:GAMMA(A) =  0',ERROR_DEF)
               RETURN
             END IF
             IF (BGAMMA == 0 ) then
               IFAIL = 3
               CALL log_mess('MATH_QTLMAP_G01EEF:GAMMA(B) =  0',ERROR_DEF)
               RETURN
             END IF
             PDF = (DGAMMA(A+B))/(AGAMMA*BGAMMA)
             PDF = PDF*(X**(A-1))*((1-X)**(B-1))
           ENDIF
#endif
           call log_mess('MATH_QTLMAP_G01EEF (OUT) ==== P:'//str(P),DEBUG_DEF)
           call log_mess('MATH_QTLMAP_G01EEF (OUT) ==== Q:'//str(Q),DEBUG_DEF)
           call log_mess('MATH_QTLMAP_G01EEF (OUT) ==== PDF:'//str(PDF),DEBUG_DEF)
           call log_mess('MATH_QTLMAP_G01EEF (OUT) ==== IFAIL:'//str(IFAIL),DEBUG_DEF)
      end subroutine MATH_QTLMAP_G01EEF


!****************************************************************************************************************

!!****f* m_qtlmap_math/MATH_QTLMAP_LOWERTAIL_BETA
!!  NAME
!!    MATH_QTLMAP_LOWERTAIL_BETA
!!  DESCRIPTION
!!   Computes the lower of the beta distribution with parameters a and b.
!!  NOTES
!!   NAG : G01EEF
!!
!!***
      subroutine MATH_QTLMAP_LOWERTAIL_BETA(X, A, B, TOL, P, IFAIL)
           double precision,intent(in)     ::  X, A, B, TOL
           double precision,intent(out)    ::  P
           integer , intent(inout)         ::  IFAIL

#ifdef HAVE_LIBNAG
           double precision                :: Q, PDF
           Q=0.d0
           PDF=0.d0

#endif
           P=0.d0
#ifdef DEBUG
           call log_mess('MATH_QTLMAP_LOWERTAIL_BETA (IN) ==== X:'//str(X),DEBUG_DEF)
           call log_mess('MATH_QTLMAP_LOWERTAIL_BETA (IN) ==== A:'//str(A),DEBUG_DEF)
           call log_mess('MATH_QTLMAP_LOWERTAIL_BETA (IN) ==== B:'//str(B),DEBUG_DEF)
#endif
           IF ( X < 0 .or. X > 1 ) then
              IFAIL = 1
              CALL log_mess('MATH_QTLMAP_LOWERTAIL_BETA:'//str(X)//' have to be between 0 and 1.',ERROR_DEF)
              RETURN
           END IF

           IF ( A < 0 )  then
             IFAIL = 2
             CALL log_mess('MATH_QTLMAP_LOWERTAIL_BETA :A is less than 0',ERROR_DEF)
             RETURN
           END IF

           IF ( B < 0 ) then
             IFAIL = 2
             CALL log_mess('MATH_QTLMAP_LOWERTAIL_BETA :B is less than 0',ERROR_DEF)
             RETURN
           END IF

#ifdef HAVE_LIBNAG
           call G01EEF(X, A, B, TOL, P, Q, PDF, IFAIL)
#else
           !The lower tail is also knwon as incomplete beta function
           P = DBETAI(X,A,B)
#endif
#ifdef DEBUG
           call log_mess('MATH_QTLMAP_LOWERTAIL_BETA (OUT) ==== P:'//str(P),DEBUG_DEF)
           call log_mess('MATH_QTLMAP_LOWERTAIL_BETA (OUT) ==== IFAIL:'//str(IFAIL),DEBUG_DEF)
#endif
      end subroutine MATH_QTLMAP_LOWERTAIL_BETA

!****************************************************************************************************************
!!****f* m_qtlmap_math/MATH_QTLMAP_G01EAF
!!  NAME
!!    MATH_QTLMAP_G01EAF
!!  DESCRIPTION
!!   Returns a one or two-tail probability for the standard Normal distribution, via the routine name.
!!
!!   The lower tail probability :
!!          P(X<=x) = 1/2*derfc(-x/sqrt(2))
!!
!!   The upper tail probability :
!!          P(X>=x) = P(X<=-x)
!!
!!   'S':P(X>=x)+P(X<=-|x|)
!!
!!   'C':P(X<=x)-P(X<=-|x|)
!!  NOTES
!!***
      function MATH_QTLMAP_G01EAF (TAIL, X, IFAIL) result (P)
           INTEGER ,intent(inout)      :: IFAIL
           double precision,intent(in) :: X
           CHARACTER(len=1),intent(in) :: TAIL

           double precision            :: P
#ifdef HAVE_LIBNAG
           double precision, external :: G01EAF
           P = G01EAF (TAIL, X, IFAIL)
#else
           double precision            :: X1,X2
           IFAIL = 0
           ! The lower tail probability : P(X<=x) = 1/2*derfc(-x/sqrt(2))
           if (TAIL == 'L' ) then
               X1=-X/SQRT_2
               P = 0.5*derfc(X1)
               return
           end if

           ! The upper tail probability : P(X>=x) = P(X<=-x)
           if ( TAIL == 'U') then
               X1=x/SQRT_2
               P = 0.5*derfc(X1)
               return
           end if

           !P(X>=x)+P(X<=-|x|)
           if ( TAIL == 'S') then
              X1=x/SQRT_2
              X2=-abs(x)/SQRT_2
              P = 0.5*derfc(X1) + 0.5*derfc(X2)
              return
           end if

           !P(X<=x)-P(X<=-|x|)
           if ( TAIL == 'C') then
             X1=-x/SQRT_2
             X2=-abs(x)/SQRT_2
             P = 0.5*derfc(X1) - 0.5*derfc(X2)
             return
           end if

           IFAIL = 1
           call stop_application('MATH_QTLMAP_G01EAF : unknown value for TAIL:'//TAIL)

#endif
           return
      end function MATH_QTLMAP_G01EAF

!****************************************************************************************************************

!!****f* m_qtlmap_math/MATH_QTLMAP_G01FAF
!!  NAME
!!    MATH_QTLMAP_G01FAF
!!  DESCRIPTION
!!   Returns the deviate associated with the given probability of the standard Normal distribution, via
!!   the routine name.
!!  NOTES
!!***
       function MATH_QTLMAP_G01FAF (TAIL, P, IFAIL) result (X)
           INTEGER ,intent(inout)      :: IFAIL
           double precision,intent(in) :: P
           CHARACTER(len=1),intent(in) :: TAIL
           double precision            :: X

#ifdef HAVE_LIBNAG
           double precision, external :: G01FAF
           x = G01FAF (TAIL, P, IFAIL)
#else
           x = dinvnorm(p)
#endif
           return

      end function MATH_QTLMAP_G01FAF


!****************************************************************************************************************

!!****f* m_qtlmap_math/MATH_QTLMAP_G01AAF
!!  NAME
!!    MATH_QTLMAP_G01AAF
!!  DESCRIPTION
!!   Calculates the mean, standard deviation, coefficients of skewness and kurtosis, and the maximum
!!   and minimum values for a set of ungrouped data. Weighting may be used.
!!  NOTES
!!***
      subroutine MATH_QTLMAP_G01AAF(N, X, IWT, WT, XMEAN, S2, S3, S4, XMIN, XMAX, WTSUM,IFAIL)
           integer , intent(in)          :: N
           integer , intent(inout)       :: IWT
           integer, intent(out)          :: IFAIL
           double precision,intent(in)   :: X(N)
           double precision,intent(inout):: WT(N)
           double precision,intent(out)  :: XMEAN, S2, S3, S4, XMIN, XMAX, WTSUM

#ifndef HAVE_LIBNAG
           integer                       :: I
           double precision              :: SOMMWTCARR,BUF,BUF2,D
#endif
#ifdef DEBUG
           call log_mess('MATH_QTLMAP_G01AAF (IN) ==== N:'//str(N),DEBUG_DEF)
          ! call log_mess('MATH_QTLMAP_G01AAF (IN) ==== X:'//str(X),DEBUG_DEF)
           call log_mess('MATH_QTLMAP_G01AAF (IN) ==== IWT:'//str(IWT),DEBUG_DEF)
#endif
#ifdef HAVE_LIBNAG
           call G01AAF(N, X, IWT, WT, XMEAN, S2, S3, S4, XMIN, XMAX, WTSUM,IFAIL)

#else

           IF (N < 1) THEN
              IFAIL = 1
              CALL log_mess('MATH_QTLMAP_G01AAF: N < 1 N:'//str(N),ERROR_DEF)
              RETURN
           END IF

           !no weight is associated with data
           IF (IWT == 0) THEN
              WT=1.d0
              IWT=1
           ELSE
              DO I=1,N
               IF ( WT(I) < 0 ) THEN
                IFAIL = 3
                CALL log_mess('MATH_QTLMAP_G01AAF: Weigth is negative for element :'//str(I),ERROR_DEF)
                RETURN
               END IF
              END DO
           END IF

           WTSUM=0.d0
           XMEAN=0.d0
           SOMMWTCARR=0.d0
           !Somme of weight and Mean
           DO I=1,N
             WTSUM = WTSUM + WT(I)
             XMEAN = XMEAN + (WT(I)*X(I))
             SOMMWTCARR = SOMMWTCARR + (WT(I)*WT(I))
           END DO

           IF ( WTSUM == 0 ) THEN
                IFAIL = 2
                CALL log_mess('MATH_QTLMAP_G01AAF: Weigth somme equal 0 ',ERROR_DEF)
                RETURN
           END IF

           XMEAN = XMEAN / WTSUM

           !Standart deviation - Coefficient of skewness - Coefficient of kurtosis
           D = WTSUM - (SOMMWTCARR / WTSUM )

           IF ( D == 0 ) THEN
                IFAIL = 2
                CALL log_mess('MATH_QTLMAP_G01AAF:'//&
                ' Can not compute Standart deviation-Coefficient of skewness- Coefficient of kurtosis',ERROR_DEF)
                RETURN
           END IF

           XMIN = X(1)
           XMAX = X(1)
           BUF= 0.d0
           S2 = 0.d0
           S3 = 0.d0
           S4 = 0.d0
           DO I=1,N
             IF ( XMIN > X(I) ) XMIN = X(I)
             IF ( XMAX < X(I) ) XMAX = X(I)
             BUF = (X(I)-XMEAN)
             BUF2 = WT(I)*BUF*BUF
             S2 = S2 + BUF2
             BUF2 = BUF2*BUF
             S3 = S3 + BUF2
             BUF2 = BUF2*BUF
             S4 = S4 + BUF2
           END DO
           S2 = S2 / D
           S2 = sqrt(S2)

           S3 = S3 / ( D*(S2**3) )
           S4 = ( S4 / ( D*(S2**4) ) ) - 3

#endif
#ifdef DEBUG
           call log_mess('MATH_QTLMAP_G01AAF (OUT) ==== XMEAN:'//str(XMEAN),DEBUG_DEF)
           call log_mess('MATH_QTLMAP_G01AAF (OUT) ==== S2:'//str(S2),DEBUG_DEF)
           call log_mess('MATH_QTLMAP_G01AAF (OUT) ==== S3:'//str(S3),DEBUG_DEF)
           call log_mess('MATH_QTLMAP_G01AAF (OUT) ==== S4:'//str(S4),DEBUG_DEF)
           call log_mess('MATH_QTLMAP_G01AAF (OUT) ==== XMIN:'//str(XMIN),DEBUG_DEF)
           call log_mess('MATH_QTLMAP_G01AAF (OUT) ==== XMAX:'//str(XMAX),DEBUG_DEF)
           call log_mess('MATH_QTLMAP_G01AAF (OUT) ==== WTSUM:'//str(WTSUM),DEBUG_DEF)
           call log_mess('MATH_QTLMAP_G01AAF (OUT) ==== IFAIL:'//str(IFAIL),DEBUG_DEF)
#endif
      end subroutine MATH_QTLMAP_G01AAF

!****************************************************************************************************************
!!****f* m_qtlmap_math/MATH_QTLMAP_G01ECF
!!  NAME
!!    MATH_QTLMAP_G01ECF
!!  DESCRIPTION
!!    Returns the lower or upper tail probability for the Chi2 distribution with real degrees of freedom, via
!!    the routine name
!!  NOTES
!!   SLATEC :
!!     * DGAMI : incomplete gamma function
!!***
      function MATH_QTLMAP_G01ECF(TAIL, X, DF, IFAIL) result(res)
         use Incomplete_Gamma
         integer , intent(out)        :: IFAIL
         real(kind=dp),intent(in)     :: X, DF
         character(len=1) ,intent(in) :: TAIL

#ifdef HAVE_LIBNAG
         real, external :: g01ecf
#else
      !   double precision,external       :: DGAMI ! SLATEC : incomplete gamma function
         real(kind=dp)                :: A,B,ans, qans,d1,d2
         integer                      :: ind
         real                         :: Xreal,res2
         real, external :: prob
#endif

         real                            :: res
#ifdef DEBUG
         call log_mess('MATH_QTLMAP_G01ECF (IN) ==== TAIL:'//TAIL,DEBUG_DEF)
         call log_mess('MATH_QTLMAP_G01ECF (IN) ==== X:'//str(X),DEBUG_DEF)
         call log_mess('MATH_QTLMAP_G01ECF (IN) ==== DF:'//str(DF),DEBUG_DEF)
#endif
         IFAIL = 0
         IF ( TAIL /= 'L' .and. TAIL /= 'U' ) THEN
             IFAIL = 1
             call log_mess('MATH_QTLMAP_G01ECF (OUT) ==== IFAIL:'//str(IFAIL),DEBUG_DEF)
             CALL log_mess('MATH_QTLMAP_G01ECF: Bad value of TAIL (U or L):'//TAIL,ERROR_DEF)
             RETURN
         END IF

         IF ( X < 0 ) THEN
             IFAIL = 2
             call log_mess('MATH_QTLMAP_G01ECF (OUT) ==== IFAIL:'//str(IFAIL),DEBUG_DEF)
             CALL log_mess('MATH_QTLMAP_G01ECF: X have to be positive or equal to 0:'//str(X),WARNING_DEF)
             RETURN
         END IF

         IF ( DF <= 0 ) THEN
             IFAIL = 3
             call log_mess('MATH_QTLMAP_G01ECF (OUT) ==== IFAIL:'//str(IFAIL),DEBUG_DEF)
             CALL log_mess('MATH_QTLMAP_G01ECF: DF have to be strictly positive:'//str(DF),WARNING_DEF)
             RETURN
         END IF

#ifdef HAVE_LIBNAG
         RES = G01ECF(TAIL, X, DF, IFAIL)
#else
        ind=int(DF)
        Xreal=X
        RES=prob(Xreal,ind)

         IF ( TAIL == 'U' ) RETURN
         IF ( TAIL == 'L') RES = 1-RES
!

#endif
#ifdef DEBUG
         call log_mess('MATH_QTLMAP_G01ECF (OUT) ==== RES:'//str(RES),DEBUG_DEF)
         call log_mess('MATH_QTLMAP_G01ECF (OUT) ==== IFAIL:'//str(IFAIL),DEBUG_DEF)
#endif
         RETURN
      end function MATH_QTLMAP_G01ECF

!****************************************************************************************************************
!!****f* m_qtlmap_math/MATH_QTLMAP_G01EBF
!!  NAME
!!    MATH_QTLMAP_G01EBF
!!  DESCRIPTION
!!   Returns the lower tail, upper tail or two-tail probability for the Student’s t-distribution with real
!!   degrees of freedom, via the routine name.
!!  NOTES
!!  HISTORY
!!   06/10/2009 : use statlib student function and stop implementation with dbetai from slatec
!!***
      function MATH_QTLMAP_G01EBF(TAIL, T, DF, IFAIL) result(res)
         integer , intent(out)        :: IFAIL
         double precision,intent(in)  :: T, DF
         character(len=1) ,intent(in) :: TAIL
#ifdef HAVE_LIBNAG
         double precision, external :: g01ebf
#endif
         double precision             :: RES

#ifdef DEBUG
         call log_mess('MATH_QTLMAP_G01EBF (IN) ==== TAIL:'//TAIL,DEBUG_DEF)
         call log_mess('MATH_QTLMAP_G01EBF (IN) ==== T:'//str(T),DEBUG_DEF)
         call log_mess('MATH_QTLMAP_G01EBF (IN) ==== DF:'//str(DF),DEBUG_DEF)
#endif
         IF ( TAIL /= 'L' .and. TAIL /= 'U' .and. TAIL /= 'C' .and. TAIL /='S' ) THEN
             IFAIL = 1
             CALL log_mess('MATH_QTLMAP_G01EBF: Bad value of TAIL (U,S,C,L):'//TAIL,ERROR_DEF)
             RETURN
         END IF

         IF ( DF < 1 ) THEN
             IFAIL = 3
             CALL log_mess('MATH_QTLMAP_G01EBF: DF have to be strictly superior to one:'//str(DF),ERROR_DEF)
             RETURN
         END IF


#ifdef HAVE_LIBNAG

         RES = G01EBF(TAIL, T, DF, IFAIL)
#else
         IF ( TAIL == 'L') RES = 1-STUDNT(T,DF,IFAIL)
         IF ( TAIL == 'U') RES = STUDNT(T,DF,IFAIL)
         IF ( TAIL == 'S' ) RES = 2*STUDNT(T,DF,IFAIL)
         IF ( TAIL == 'C' ) RES = 1 - 2*STUDNT(T,DF,IFAIL)

#endif
#ifdef DEBUG
        call log_mess('MATH_QTLMAP_G01EBF (OUT) ==== RES:'//str(RES),DEBUG_DEF)
        call log_mess('MATH_QTLMAP_G01EBF (OUT) ==== IFAIL:'//str(IFAIL),DEBUG_DEF)
#endif
      end function MATH_QTLMAP_G01EBF


!****************************************************************************************************************
!!****f* m_qtlmap_math/MATH_QTLMAP_G03ACF
!!  NAME
!!    MATH_QTLMAP_G05EHF
!!  DESCRIPTION
!!   Performs a canonical variate (canonical discrimination) analysis.
!!  INPUT
!!
!!     WEIGHT    : 'U' no weights ar used . ’W’ or ’V’, weights are used and must be supplied in WT.
!!                 Constraint: WEIGHT 1⁄4 ’U’; ’W’ or ’V’.
!!                 In the case of WEIGHT 1⁄4 ’W’, the weights are treated as frequencies and the effective number of
!!                 observations is the sum of the weights. If WEIGHT 1⁄4 ’V’, the weights are treated as being
!!                 inversely proportional to the variance of the observations and the effective number of observations is
!!                 the number of observations with non-zero weights.
!!
!!     N         : the number of observations. Constraint: N >= NX + NG.
!!     M         : the total number of variables. Constraint: M >= NX.
!!     X(LDX,M)  : X(i,j) must contain the ith observation for the jth variable, 1<=i<=N, 1<=j<=M
!!     LDX       : the first dimension of the array X LDX >= N
!!     ISX(M)    : indicates whether or not the jth variable is to be included in the analysis.
!!                 If ISX(j) > 0, then the variables contained in the jth column of X is included in the canonical variate analysis,
!!     NX        : the number of variables in the analysis, nx .Constraint : NX >= 1
!!     ING(N)    : ING(i) indicates which group the ith observation is in, for i=1,2,...,n. The effective
!!                 number of groups is the number of groups with non-zero membership.
!!                 Constraint 1<=NG(i)<= NG i=1,..,N
!!     NG        : the number of groups, ng . Constraint: NG >=2
!!     WT        : size(WT)>=N if WEIGHT='W' or 'V' , 1 otherwise.
!!                 if WEIGHT=’W’ or ’V’, then the first n elements of WT must contain the weights to used in the analysis.
!!                 if W(i)=0.0, then the ith observation is not included in the analysis.
!!                 if WEIGHT=’U’, then WT is not referenced.
!!                 Constraints:
!!                             WT(i)>=0.0,for i=1,2,..,N   ,
!!                             SOMME(WT(i))>=NX + effective number of groups , i=1,...,N
!!     LDCVM     : the dimension of the array CVM as declared in the program. Constraint LDCVM >=NG
!!     LDE       : the first dimension of the array E. Constraint: LDE >= min(NX,NG-1).
!!     LDCVX     : the first dimension of the array CVX
!!     TOL       : the value of TOL is used to decide if the variables are of full rank and, if not, what is the
!!                 rank of the variables. The smaller the value of TOL the stricter the criterion for selecting the
!!                 singular value decomposition. If a non-negative value of TOL less than machine precision is
!!                 entered, then the square root of machine precision is used instead.
!!
!!  OUTPUT
!!
!!    NIG(NG)       : NIG(j) gives the number of observations in group j, for j=1,2,...,ng.
!!    CVM(LDCVM,NX) : CVM(i,j) contains the mean of the jth canonical variate for the ith group, for
!!                   i=1,2,...,ng ; j=1,2,...,l the remaining columns, if any, are used as workspace.
!!    E(LDE,6)      : the statistics of the canonical variate analysis.
!!                   E(i,1), the canonical correlations, deltai , for i=1,2,...,l.
!!                   E(i,2), the eigenvalues of the within-group sum of squares matrix ,lambda2, for i=1,...,l.
!!                   E(i,3), the proportion of variation explained by the ith canonical variate, for i=1,...,l.
!!                   E(i,4), the Chi2 statistic for the ith canonical variate, for i=1,...,l.
!!                   E(i,5), the degrees of freedom for Chi2 statistic for the ith canonical variate, for i=1,...,l.
!!                   E(i,6), he significance level for the Chi2 statistic for the ith canonical variate, for i=1,...,l.
!!    NCV           :  he number of canonical variates, l. This will be the minimum of NG-1 and the rank of X.
!!    CVX(LDCVX,NG-1) : the canonical variate loadings. CVX(i,j) contains the loading coefficient for the ith
!!                      variable on the jth canonical variate, for i=1,...,NX ; j=1,...,l; the remaining columns, if
!!                      any, are used as workspace.
!!    IRANKX        :  the rank of the dependent variables.
!!                     If the variables are of full rank then IRANKX 1⁄4 NX.
!!                     If the variables are not of full rank then IRANKX is an estimate of the rank of the dependent
!!                    variables. IRANKX is calculated as the number of singular values greater than TOLÂ (largest singular value).
!!
!!
!!  NOTES
!!***
    subroutine  MATH_QTLMAP_G03ACF(WEIGHT, N, M, X, LDX, ISX, NX, ING, NG, WT, NIG, CVM, &
                LDCVM, E, LDE, NCV, CVX, LDCVX, TOL, IRANKX,IFAIL)
        CHARACTER(LEN=1)                   ,INTENT(in)     :: WEIGHT
        INTEGER                            ,INTENT(in)     :: N, M, LDX, NX, NG
        INTEGER                            ,INTENT(in)     :: LDCVM, LDE, LDCVX
        DOUBLE PRECISION  ,DIMENSION(LDX,M),INTENT(in)     :: X
        INTEGER           ,DIMENSION(M)    ,INTENT(in)     :: ISX
        INTEGER           ,DIMENSION(N)    ,INTENT(in)     :: ING
        DOUBLE PRECISION  ,DIMENSION(:)    ,INTENT(in)     :: WT
        DOUBLE PRECISION                   ,INTENT(in)     :: TOL
   !output
        INTEGER           ,DIMENSION(NG)   ,INTENT(out)    :: NIG
        DOUBLE PRECISION  ,DIMENSION(LDCVM,NX),INTENT(out) :: CVM
        DOUBLE PRECISION  ,DIMENSION(LDE,6),INTENT(out)    :: E
        INTEGER                            ,INTENT(out)    :: NCV
        DOUBLE PRECISION  ,DIMENSION(LDCVX,NG-1),INTENT(out):: CVX

        INTEGER                            ,INTENT(out)    :: IRANKX,  IFAIL



#ifdef HAVE_LIBNAG
      INTEGER                          :: IWK,STAT
      DOUBLE PRECISION  ,DIMENSION(:),ALLOCATABLE  :: WK     ! workspace
      external :: G03ACF

      IF (NX >= NG -1) THEN
         IWK = 1 + N*NX + max ( (5*(NX-1)) + ((NX+1)*NX) ,N)
      ELSE
         IWK = 1 + N*NX + max ( (5*(NX-1)) + ((NG-1)*NX),N)
      ENDIF

      allocate(WK(IWK),STAT=stat)

      if (stat /= 0 ) then
        call log_mess("Can not allocate Workspace for G03ACF",ERROR_DEF)
        stop
      end if

      call G03ACF(WEIGHT, N, M, X, LDX, ISX, NX, ING, NG, WT, NIG, CVM, &
                LDCVM, E, LDE, NCV, CVX, LDCVX, TOL, IRANKX, WK, IWK, &
                IFAIL)
      deallocate(WK);
#else
        INTEGER  :: I,J,IPRINT,IERR,IG,M2,K,MATZ,IP
      INTEGER  :: JP,MPRIME,LDVL,LDVR,LWORK,INFO,NPRIME,CORR(M),BIG
      INTEGER, DIMENSION(NG) :: NOG
      DOUBLE PRECISION  :: MEAN_VAR(M),MEAN_GROUP(NG,M)
      DOUBLE PRECISION  :: WT2(N),WT22(N),WR2(M),NIGR(NG)
      DOUBLE PRECISION  ,DIMENSION(:,:)  ,ALLOCATABLE  :: X2,WITHIN,BETWEEN
      INTEGER           ,DIMENSION(N)        :: ING2

      DOUBLE PRECISION  :: WR(M),WI(M),A(M,M),Z(M,M),VL(M,M),VR(M,M),WORK(100*M)
      CHARACTER :: JOBVL,JOBVR
      EXTERNAL :: DGEEV


         IF ( NX < 1 ) THEN
         print *,"MATH_QTLMAP_G03ACF : NX < 1"
         !call log_mess("MATH_QTLMAP_G03ACF : NX < 1",ERROR_DEF);
         IFAIL=1
         RETURN
      END IF

      IF ( NG < 2 ) THEN
         print *,"MATH_QTLMAP_G03ACF : NG < 2"
         !call log_mess("MATH_QTLMAP_G03ACF : NG < 2",ERROR_DEF);
         IFAIL=1
         RETURN
      END IF

      IF ( M < NX ) THEN
         print *,"MATH_QTLMAP_G03ACF : M < NX"
         !call log_mess("MATH_QTLMAP_G03ACF : M < NX",ERROR_DEF);
         IFAIL=1
         RETURN
      END IF

     IF ( N < (NX+NG) ) THEN
         print *,"MATH_QTLMAP_G03ACF : N < NX+NG"
         !call log_mess("MATH_QTLMAP_G03ACF : N < NX+NG",ERROR_DEF);
         IFAIL=1
         RETURN
      END IF

     IF ( LDX < N ) THEN
         print *,"MATH_QTLMAP_G03ACF : LDX < N"
         !call log_mess("MATH_QTLMAP_G03ACF : LDX < N",ERROR_DEF);
         IFAIL=1
         RETURN
      END IF

     IF ( LDCVX < NX ) THEN
         print *,"MATH_QTLMAP_G03ACF : LDCVX < NX"
         !call log_mess("MATH_QTLMAP_G03ACF : LDCVX < NX",ERROR_DEF);
         IFAIL=1
         RETURN
      END IF

     IF ( LDCVM < NG ) THEN
         print *,"MATH_QTLMAP_G03ACF : LDCVM < NG"
         !call log_mess("MATH_QTLMAP_G03ACF : LDCVM < NG",ERROR_DEF);
         IFAIL=1
         RETURN
      END IF

     IF ( LDE < min(NX,NG-1) ) THEN
         print *,"MATH_QTLMAP_G03ACF : LDE < min(NX,NG-1)"
         !call log_mess("MATH_QTLMAP_G03ACF : LDE < min(NX,NG-1)",ERROR_DEF);
         IFAIL=1
         RETURN
      END IF

     IF ( WEIGHT/='U' .and. WEIGHT/='W' .and. WEIGHT/='V' ) THEN
         print *,"MATH_QTLMAP_G03ACF : WEIGHT /= 'U' or 'V' or 'W' "
         !call log_mess("MATH_QTLMAP_G03ACF : WEIGHT /= 'U' or 'V' or 'W' ",ERROR_DEF);
         IFAIL=1
         RETURN
      END IF

      IF ( TOL < 0.d0 ) THEN
         print *,"MATH_QTLMAP_G03ACF : TOL < 0.0"
         !call log_mess("MATH_QTLMAP_G03ACF : TOL < 0.0",ERROR_DEF);
         IFAIL=1
         RETURN
      END IF
    ! IFAIL = 2
       IF ( ( WEIGHT=='V' .or. WEIGHT=='W') .and. minval(WT)<0.d0  ) THEN
         print *,"MATH_QTLMAP_G03ACF : WEIGHT:"//WEIGHT//" AND ONE VALUE OF WT < 0.0"
         !call log_mess("MATH_QTLMAP_G03ACF : WEIGHT:"//WEIGHT//" AND ONE VALUE OF WT < 0.0",ERROR_DEF);
         IFAIL=2
         RETURN
      END IF

    IFAIL=0


    !Weights can be used with the analysis, in which case the weighted means are subtracted from each column
    !and then each row is scaled by an amount square(wi) , where wi is the weight for the ith observation (row).

     NPRIME=0

!     DO I=1,N
!        X2(I,:) = X(I,:)*wt(i)
!        IF (WT(I)>0) NPRIME=NPRIME+WT(I)
!     END DO

     IF ( WEIGHT == 'U') THEN
       WT2 = 1.d0
     ELSE IF ( WEIGHT == 'W' ) THEN
       WT2 = WT
     ELSE IF ( WEIGHT == 'V' ) THEN
       WT2 = 0.d0
       print *," ** ko ** WEIGHT=V"
       stop
     END IF
     IPRINT=3
     IERR=0

     NOG=0
     !*********** Calcul du nombre d observation dans chaque groupe
     NIG=0
     NIGR=0.d0
     !
     DO I=1,N
       IF ( ING(I) > 0 .AND. (WT(I)>0) ) THEN
         NIG(ING(I))= NIG(ING(I)) + 1
         NIGR(ING(I))= NIGR(ING(I)) + WT2(I)
       END IF
     END DO

     NPRIME=N

     !*********** calcul des moyennes de chaque variable
     MEAN_VAR=0d0
     JP=0
     DO J=1,M
      IF (ISX(J)>0) THEN
       JP=JP+1
       DO I=1,N
        IF ( ING(I) > 0 ) THEN
          MEAN_VAR(JP) = X(I,JP)*WT2(I) + MEAN_VAR(JP)
        END IF
       END DO
       MEAN_VAR(JP) = MEAN_VAR(JP) / sum(WT2(:N))
      END IF
     END DO

     !NB VARIABLES
     MPRIME=JP

     ALLOCATE (X2(N,MPRIME))
     ALLOCATE (WITHIN(MPRIME,MPRIME))
     ALLOCATE (BETWEEN(MPRIME,MPRIME))


     JP=0
     DO J=1,M
      IF (ISX(J)>0) THEN
       JP=JP+1
       IP=0
       DO I=1,N
        IF ( ING(I) > 0 ) THEN
          IP=IP+1
          X2(IP,JP) = X(I,J)
          ING2(IP)= ING(I)
          WT22(IP) = WT2(I)
        END IF
       END DO
      END IF
     END DO

     IP=-1
     JP=-1

     !*********** calcul des moyennes de chaque groupes et variables
     MEAN_GROUP=0d0
     DO I=1,NPRIME
      DO J=1,MPRIME
         IF(ING2(I) == 0) PRINT *,I,ING2(I),X(I,J),X2(I,J),ING(I)
         MEAN_GROUP(ING2(I),J)= X2(I,J)*WT22(I) + MEAN_GROUP(ING2(I),J)
      END DO
     END DO

     DO J=1,MPRIME
       DO I=1,NG
        MEAN_GROUP(I,J)= MEAN_GROUP(I,J) / NIGR(I)
       END DO
     END dO
!
     WITHIN=0.d0
     DO I=1,MPRIME
      DO J=1,MPRIME
        DO IG=1,NG
         DO K=1,NPRIME
          IF ( ING2(K) == IG ) THEN
            WITHIN(I,J)= WITHIN(I,J) + (X2(K,I) - MEAN_GROUP(ING2(K),I))*(X2(K,J) - MEAN_GROUP(ING2(K),J))*WT22(K)
          END IF
         END DO
        END DO
      END DO
     END DO

     !*********** Calcul de Between k (k=1,..,NG) : matrice de variance covariance
     BETWEEN=0.d0

     DO I=1,MPRIME
       DO J=1,MPRIME
         DO K=1,NG
          BETWEEN(I,J)= BETWEEN(I,J) + ( NIGR(K) )* &
          (MEAN_GROUP(K,I) - MEAN_VAR(I)) * (MEAN_GROUP(K,J) - MEAN_VAR(J))
         END DO
       END DO
     END DO

     call ginv1(WITHIN,M,MPRIME,TOL,INFO)
     A = matmul(WITHIN,BETWEEN)

     WR=0
     Wi=0

     JOBVL='N'
     JOBVR='V'
     LDVL=M
     LDVR=M
     LWORK=10*N
     CALL DGEEV(JOBVL,JOBVR,M,A,M,WR,WI,VL,LDVL,VR,LDVR,WORK,LWORK,INFO)
     IF ( INFO /= 0 ) THEN
       print *,"DGEEV fail in G03ACF"
     END IF

     NCV=0 ! nombre de valeur propre > tol
     WR2=WR
     DO I=1,M
       BIG=0
       DO J=1,M
        IF ( WR2(J) > TOL ) THEN
         IF (BIG > 0) THEN
          IF (WR2(BIG) < WR2(J)) BIG = J
         ELSE
            BIG = J
         END IF
        END IF
       END DO
       IF (BIG > 0 ) THEN
         WR2(BIG) = 0.d0
         NCV=NCV+1
         CORR(NCV)=BIG
       ELSE
        EXIT
       END IF
     END DO

     NCV = min(NCV,NG-1)

     !Canonical correlation
      DO I=1,NCV
       E(I,1) = sqrt(WR(CORR(I)) / (1.d0 + WR(CORR(I))))
       E(I,2) = WR(CORR(I))
       CVX(:M,I) = VR(:M,CORR(I))
     END DO

     DEALLOCATE (X2)
     DEALLOCATE (WITHIN)
     DEALLOCATE (BETWEEN)

#endif
    end subroutine MATH_QTLMAP_G03ACF


   subroutine test_G03ACF
     INTEGER :: N=9
     INTEGER :: M=3
     INTEGER :: NX=3
     INTEGER :: NG=3
     CHARACTER(LEN=1) :: WEIGHT='U'

     double precision    TOL
     INTEGER IFAIL, IRX, J, NCV
    double precision  ::CVM(9,9), CVX(9,9), E(9,6),WK(50)
    INTEGER ISX(6), NIG(9)


     DOUBLE PRECISION :: X(9,3),WT(9)
     INTEGER          :: ING(9),I

     data (X(1,i),i=1,3) /13.3,10.6,21.2/
     data (X(2,i),i=1,3) /13.6,10.2,21.0/
     data (X(3,i),i=1,3) /14.2,10.7,21.1/
     data (X(4,i),i=1,3) /13.4,9.4,21.0/
     data (X(5,i),i=1,3) /13.2,9.6,20.1/
     data (X(6,i),i=1,3) /13.9,10.4,19.8/
     data (X(7,i),i=1,3) /12.9,10.0,20.5/
     data (X(8,i),i=1,3) /12.2,9.9,20.7/
     data (X(9,i),i=1,3) /13.9,11.0,19.1/
     data ING /1,2,3,1,2,3,1,2,3/

     ISX=1
     TOL=1.d-5
     IFAIL=0

     PRINT *,X(1,:)
     PRINT *,ING

     CALL MATH_QTLMAP_G03ACF(WEIGHT,N,M,X,9,ISX,NX,ING,NG,WT,&
     NIG,CVM,3,E,3,NCV,CVX,3,TOL,IRX,IFAIL)

    WRITE (*,*)
    WRITE (*,99999) 'Rank of X = ', IRX
    WRITE (*,*)
           WRITE (*,*) 'Canonical     Eigenvalues Percentage    CHISQ       DF SIG'
           WRITE (*,*) 'Correlations              Variation'
           DO 60 I = 1, NCV
              WRITE (*,99998) (E(I,J),J=1,6)
   60      CONTINUE
           WRITE (*,*)
           WRITE (*,*) 'Canonical Coefficients for X'
           DO 80 I = 1, NX
              WRITE (*,99997) (CVX(I,J),J=1,NCV)
   80      CONTINUE
           WRITE (*,*)
           WRITE (*,*) 'Canonical variate means'
           DO 100 I = 1, NG
              WRITE (*,99997) (CVM(I,J),J=1,NCV)
  100      CONTINUE

       STOP

99999 FORMAT (1X,A,I2)
99998 FORMAT (1X,2F12.4,F11.4,F10.4,F8.1,F8.4)
99997 FORMAT (1X,5F9.4)

   end subroutine test_G03ACF


!****************************************************************************************************************
!!****f* m_qtlmap_math/MATH_QTLMAP_G05EHF
!!  NAME
!!    MATH_QTLMAP_G05EHF
!!  DESCRIPTION
!!   Performs a pseudo-random permutation of a vector of integers.
!!  NOTES
!!   Use order-pack module m_ctrper :
!!   Permute an array randomly, but leaving elements close
!!   to their initial locations (nearbyness is controled by PCLS).
!!***
    subroutine MATH_QTLMAP_G05EHF(INDEX, N, IFAIL)
#ifndef HAVE_LIBNAG
      use m_ctrper
#endif
      INTEGER  , INTENT(IN)                 :: N
      INTEGER  , DIMENSION(N),INTENT(INOUT) :: INDEX
      INTEGER  , INTENT(INOUT)              :: IFAIL

#ifdef HAVE_LIBNAG
      external :: G05EHF
      call G05EHF(INDEX, N, IFAIL)
#else
      integer                               :: i
      IF ( N < 1 ) THEN
        call log_mess("MATH_QTLMAP_G05EHF Dimension ARRAY have to be greater or equal than 1.",ERROR_DEF)
        IFAIL=1
        RETURN
      END IF

      call ctrper(INDEX,1.0)
#endif

    end subroutine MATH_QTLMAP_G05EHF

!****************************************************************************************************************
!!****f* m_qtlmap_math/MATH_QTLMAP_M01CAF
!!  NAME
!!    MATH_QTLMAP_M01CAF
!!  DESCRIPTION
!!   Rearranges a vector of real numbers into ascending or descending order.
!!  NOTES
!!  ORDERPACK
!!***
      subroutine MATH_QTLMAP_M01CAF(RV, M1, M2, ORDER, IFAIL)
#ifndef HAVE_LIBNAG
                use m_refsor
#endif
                integer , intent(in)         :: M1, M2
                integer , intent(out)        :: IFAIL
                double precision ,intent(inout)  :: RV(M2)
                character(len=1) ,intent(in) :: ORDER
#ifdef HAVE_LIBNAG
#else
                double precision , dimension(M2-M1+1):: DX,DY
                integer i
#endif

		          IF ( M1 < 1 ) THEN
		             IFAIL = 1
		             CALL log_mess('MATH_QTLMAP_M01CAF: M1 < 1 :'//str(M1),ERROR_DEF)
		             RETURN
		          END IF

		          IF ( M2 < 1 ) THEN
		             IFAIL = 1
		             CALL log_mess('MATH_QTLMAP_M01CAF: M2 < 1 :'//str(M2),ERROR_DEF)
		             RETURN
		          END IF

		          IF ( M1 > M2 ) THEN
		             IFAIL = 1
		             CALL log_mess('MATH_QTLMAP_M01CAF: M1 > M2  M1:'//str(M1)//' M2:'//str(M2),ERROR_DEF)
		             RETURN
		          END IF

		          IF ( ORDER /= 'A' .and. ORDER /= 'D' ) THEN
		             IFAIL = 2
		             CALL log_mess('MATH_QTLMAP_M01CAF: Bad value of ORDER (A,D):'//ORDER,ERROR_DEF)
		             RETURN
		         END IF

#ifdef HAVE_LIBNAG
                 call M01CAF(RV, M1, M2, ORDER, IFAIL)
#else
                 call refsor(RV(M1:M2))

                 IF ( ORDER == 'D' .or. ORDER == 'd' ) THEN
                     RV(M1:M2) = RV(M2:M1:-1)
                 END IF

#endif


      end subroutine MATH_QTLMAP_M01CAF

!****************************************************************************************************************
!!****f* m_qtlmap_math/MATH_QTLMAP_M01DAF
!!  NAME
!!    MATH_QTLMAP_M01DAF
!!  DESCRIPTION
!!    Ranks a vector of real numbers in ascending or descending order.
!!  NOTES
!!   ORDERPACK
!!***
      subroutine MATH_QTLMAP_M01DAF(RV, M1, M2, ORDER, IRANK, IFAIL)
#ifndef HAVE_LIBNAG
          use m_mrgrnk
#endif
          integer , intent(in)             :: M1, M2
          integer , intent(out)            :: IRANK(M2)
          integer , intent(inout)          :: IFAIL
          double precision ,intent(in)     :: RV(M2)
          character(len=1) ,intent(in)     :: ORDER

#ifdef HAVE_LIBNAG
          external :: M01DAF
#ifdef DEBUG_PRINT
          PRINT * ,' --->  M01DAF ** NAG **'
#endif
          call M01DAF(RV, M1, M2, ORDER, IRANK, IFAIL)
#else
          double precision ,dimension(1:(M2-M1+1)) :: XDONT
          integer,dimension(M2) :: RANK
          integer :: i
#ifdef DEBUG_PRINT
          PRINT * ,' --->  IM01DAF'
#endif
          IFAIL = 0
          !On entry: elements M1 to M2 of RV must contain real values to be ranked.

          !Constraintes M1 > 0 and M2 >=M1
          if ( M1 <= 0 .or. M2 < M1 ) then
            IFAIL=1
            return
          end if
          !Constraintes ORDER='A' or 'D'
          if ( ORDER /= 'A' .and. ORDER /='D') then
            IFAIL=2
            return
          end if
         ! On exit: elements M1 to M2 of IRANK contain the ranks of the corresponding elements of RV.
         ! Note that the ranks are in the range M1 to M2: thus, if RV i is the first element in the rank order,
         ! IRANK i is set to M1.

         XDONT = RV(M1:M2)

         call mrgrnk (XDONT, RANK)

         if ( ORDER == 'D') then
            RANK(:) = RANK(M2-M1+1:1:-1)
         end if
         !Rank contient l'indice du rang d'un element du tableau
         do i=1,M2-M1+1
           IRANK(RANK(i)+M1-1)=i
         end do

#endif
      end subroutine MATH_QTLMAP_M01DAF


!****************************************************************************************************************
!!  NAME
!!    IM01DBF
!!  DESCRIPTION
!!    ranks a vector of integer numbers in ascending or descending order.
!!  NOTES
!!***
      subroutine MATH_QTLMAP_M01DBF (IV, M1, M2, ORDER, IRANK, IFAIL)
#ifndef HAVE_LIBNAG
          use m_mrgrnk
#endif
          integer , intent(in)             :: M1, M2
          integer , intent(out)            :: IRANK(M2)
          integer , intent(inout)          :: IFAIL
          integer , intent(in)             :: IV(M2)
          character(len=1) ,intent(in)     :: ORDER
#ifdef HAVE_LIBNAG
          external :: M01DBF
#ifdef DEBUG_PRINT
          PRINT * ,' --->  IM01DBF ** NAG **'
#endif
          call M01DBF(IV, M1, M2, ORDER, IRANK, IFAIL)
#else
          integer ,dimension(1:(M2-M1+1)) :: XDONT
          integer,dimension(M2) :: RANK
          integer :: i
#ifdef DEBUG_PRINT
          PRINT * ,' --->  IM01DBF'
#endif
          IFAIL = 0
          !On entry: elements M1 to M2 of RV must contain real values to be ranked.

          !Constraintes M1 > 0 and M2 >=M1
          if ( M1 <= 0 .or. M2 < M1 ) then
            IFAIL=1
            return
          end if
          !Constraintes ORDER='A' or 'D'
          if ( ORDER /= 'A' .and. ORDER /='D') then
            IFAIL=2
            return
          end if
         ! On exit: elements M1 to M2 of IRANK contain the ranks of the corresponding elements of RV.
         ! Note that the ranks are in the range M1 to M2: thus, if RV i is the first element in the rank order,
         ! IRANK i is set to M1.

         XDONT = IV(M1:M2)

         call mrgrnk (XDONT, RANK)

         if ( ORDER == 'D') then
            RANK(:) = RANK(M2-M1+1:1:-1)
         end if
         IRANK=0

         !Rank contient l'indice du rang d'un element du tableau
         do i=1,M2-M1+1
          IRANK(RANK(i)+M1-1)=i
         end do

#endif
      end subroutine MATH_QTLMAP_M01DBF


Subroutine char_mrgrnk (XDONT, IRNGT)
! __________________________________________________________
!   MRGRNK = Merge-sort ranking of an array
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
! __________________________________________________________
! __________________________________________________________
      character(len=*), Dimension (:), Intent (In)  :: XDONT
      Integer, Dimension (:), Intent (Out) :: IRNGT
! __________________________________________________________
      character(len=256) :: XVALA, XVALB
!
      Integer, Dimension (SIZE(IRNGT)) :: JWRKT
      Integer :: LMTNA, LMTNC, IRNG1, IRNG2
      Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
!
      NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
      Select Case (NVAL)
      Case (:0)
         Return
      Case (1)
         IRNGT (1) = 1
         Return
      Case Default
         Continue
      End Select
!
!  Fill-in the index array, creating ordered couples
!
      Do IIND = 2, NVAL, 2
         If (lle(XDONT(IIND-1),XDONT(IIND))) Then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         Else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         End If
      End Do
      If (Modulo(NVAL, 2) /= 0) Then
         IRNGT (NVAL) = NVAL
      End If
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      LMTNA = 2
      LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      Do
         If (NVAL <= 2) Exit
!
!   Loop on merges of A and B into C
!
         Do IWRKD = 0, NVAL - 1, 4
            If ((IWRKD+4) > NVAL) Then
               If ((IWRKD+2) >= NVAL) Exit
!
!   1 2 3
!
               If (lle(XDONT(IRNGT(IWRKD+2)),XDONT(IRNGT(IWRKD+3)))) Exit
!
!   1 3 2
!
               If (lle(XDONT(IRNGT(IWRKD+1)),XDONT(IRNGT(IWRKD+3)))) Then
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
!
!   3 1 2
!
               Else
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               End If
               Exit
            End If
!
!   1 2 3 4
!
            If (lle(XDONT(IRNGT(IWRKD+2)),XDONT(IRNGT(IWRKD+3)))) Cycle
!
!   1 3 x x
!
            If (lle(XDONT(IRNGT(IWRKD+1)),XDONT(IRNGT(IWRKD+3)))) Then
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               If (lle(XDONT(IRNG2),XDONT(IRNGT(IWRKD+4)))) Then
!   1 3 2 4
                  IRNGT (IWRKD+3) = IRNG2
               Else
!   1 3 4 2
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               End If
!
!   3 x x x
!
            Else
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               If (lle(XDONT(IRNG1),XDONT(IRNGT(IWRKD+4)))) Then
                  IRNGT (IWRKD+2) = IRNG1
                  If (lle(XDONT(IRNG2),XDONT(IRNGT(IWRKD+4)))) Then
!   3 1 2 4
                     IRNGT (IWRKD+3) = IRNG2
                  Else
!   3 1 4 2
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  End If
               Else
!   3 4 1 2
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               End If
            End If
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 4
         Exit
      End Do
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is doubled.
!
      Do
         If (LMTNA >= NVAL) Exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
!
!   Loop on merges of A and B into C
!
         Do
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JINDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            If (IWRKF >= NVAL) Then
               If (JINDA >= NVAL) Exit
               IWRKF = NVAL
            End If
            IINDA = 1
            IINDB = JINDA + 1
!
!   Shortcut for the case when the max of A is smaller
!   than the min of B. This line may be activated when the
!   initial set is already close to sorted.
!
!          IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
!
!  One steps in the C subset, that we build in the final rank array
!
!  Make a copy of the rank array for the merge iteration
!
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
!
            XVALA = XDONT (JWRKT(IINDA))
            XVALB = XDONT (IRNGT(IINDB))
!
            Do
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               If (XVALA > XVALB) Then
                  IRNGT (IWRK) = IRNGT (IINDB)
                  IINDB = IINDB + 1
                  If (IINDB > IWRKF) Then
!  Only A still with unprocessed values
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                     Exit
                  End If
                  XVALB = XDONT (IRNGT(IINDB))
               Else
                  IRNGT (IWRK) = JWRKT (IINDA)
                  IINDA = IINDA + 1
                  If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                  XVALA = XDONT (JWRKT(IINDA))
               End If
!
            End Do
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      End Do
!
      Return
!
End Subroutine char_mrgrnk

!****************************************************************************************************************
!!  NAME
!!    IM01DCF
!!  DESCRIPTION
!!    Ranks a vector of character data in ASCII or reverse ASCII order of a specified substring.
!!
!!******************************************
      subroutine MATH_QTLMAP_M01DCF(CH, M1, M2, L1, L2, ORDER, IRANK, IFAIL)
                   INTEGER        M1, M2, L1, L2, IRANK(M2), IFAIL
                   CHARACTER*(*)  CH(M2)
                   CHARACTER*1    ORDER
#ifdef HAVE_LIBNAG
          external :: M01DCF
#ifdef DEBUG_PRINT
          PRINT * ,' --->  IM01DBF ** NAG **'
#endif
          call M01DCF(CH, M1, M2, L1, L2, ORDER, IRANK, IFAIL)
#else
          CHARACTER(LEN=256) ,dimension(1:(M2-M1+1)) :: XDONT
          integer,dimension(M2) :: RANK
          integer :: i
           IFAIL = 0
          !On entry: elements M1 to M2 of RV must contain real values to be ranked.

          !Constraintes M1 > 0 and M2 >=M1
          if ( M1 <= 0 .or. M2 < M1 ) then
            IFAIL=1
            return
          end if

          IF ( (L2<1) .or. (L1<1) .or. (L1>L2) ) THEN
            IFAIL=1
            RETURN
          END IF

          IF ( SIZE(CH) >= 1) THEN
              IF ( L2 > LEN(CH(1))) THEN
                IFAIL=1
                RETURN
              END IF
          END IF
          !Constraintes ORDER='A' or 'D'
          if ( ORDER /= 'A' .and. ORDER /='R') then
            IFAIL=2
            return
          end if
          DO I=M1,M2
              IF ( len(CH(I)) > 255 ) THEN
                IFAIL=3
                RETURN
              END IF
          END DO

         ! On exit: elements M1 to M2 of IRANK contain the ranks of the corresponding elements of RV.
         ! Note that the ranks are in the range M1 to M2: thus, if RV i is the first element in the rank order,
         ! IRANK i is set to M1.

         DO I=M1,M2
             XDONT(I-M1+1)=CH(I)(L1:L2)
         END DO

         call char_mrgrnk (XDONT, RANK)

         if ( ORDER == 'R') then
            RANK(:) = RANK(M2-M1+1:1:-1)
         end if
         IRANK=0

         !Rank contient l'indice du rang d'un element du tableau
         do i=1,M2-M1+1
          IRANK(RANK(i)+M1-1)=i
         end do

#endif

      end subroutine MATH_QTLMAP_M01DCF


!! Sample
!! Data Ranks
!!  34    9
!!  44    8
!!  89    2
!!  64    7
!!  69    4
!!  69    5
!!  23   10
!!   1   12
!! 999    1
!!  65    6
!!  22   11
!!  76    3

!****************************************************************************************************************
!!  NAME
!!    IM01ZAF
!!  DESCRIPTION
!!    Inverts a permutation, and hence converts a rank vector to an index vector, or vice versa.
!!    There are two common ways of describing a permutation using an integer vector IPERM.
!!    The first uses ranks: IPERMi holds the position to which the ith data element should be moved in order to sort the data; in other words its rank in the sorted order.
!!    The second uses indices: IPERMi holds the current position of the data element which would occur in ith position in sorted order.
!!    For example, given the values
!!     3.5 5.9 2.9 0.5
!!
!!     to be sorted in ascending order, the ranks would be
!!     3 4 2 1
!!
!!     and the indices would be
!!     4 3 1 2
!!
!!  PARAMETERS
!!
!! IPERM(M2) – INTEGER arrayInput/Output
!!    On entry: elements M1 to M2 of IPERM must contain a permutation of the integers M1 to M2.
!!
!!    On exit: these elements contain the inverse permutation of the integers M1 to M2.
!!
!!    M1 and M2 must specify the range of elements used in the array IPERM and the range of values in the permutation, as specified under IPERM.
!!
!!    Constraint:
!!      0<M1≤M2.
!!
!!***
      subroutine MATH_QTLMAP_M01ZAF ( IPERM, M1, M2, IFAIL)
         INTEGER                 , INTENT(IN)    :: M1, M2
         INTEGER  ,DIMENSION(M2) , INTENT(INOUT) :: IPERM
         INTEGER                   , INTENT(OUT) :: IFAIL
#ifdef HAVE_LIBNAG
          external :: M01ZAF
          INTEGER                                 :: I
#ifdef DEBUG_PRINT
          PRINT * ,' --->  M01ZAF ** NAG **'
#endif

          call M01ZAF ( IPERM, M1, M2, IFAIL)
#ifdef DEBUG_PRINT

          PRINT *,"************  ENV NAG *****************"
          DO I=1,M2
            PRINT *,'IPERM(',I,')=',IPERM(I)
          END DO
#endif
#else
         LOGICAL ,DIMENSION(M2)  :: PRES
         INTEGER ,DIMENSION(M2)  :: INDEX_RANK
         INTEGER                 :: I,J
#ifdef DEBUG_PRINT
          PRINT * ,' --->  IM01ZAF'
#endif

         IF ( (M2<1) .or. (M1<1) .or. (M1>M2) ) THEN
           IFAIL=1
           RETURN
         END IF

         PRES=.false.

         DO I=M1,M2
           IF ( (IPERM(I)>M2) .or. IPERM(I)<M1  ) THEN
             IFAIL = 2
             RETURN
           END IF
           IF ( PRES(IPERM(I)) ) THEN
             IFAIL = 3
             RETURN
           END IF
           PRES(IPERM(I)) = .true.
         END DO

         INDEX_RANK=IPERM
         DO I=M1,M2 ! for each rank
          ! find the index of IPERM that contains the minimum
           DO J=M1,M2
            IF ( I == IPERM(J) ) THEN
                INDEX_RANK(I) = J
                EXIT
            END IF
           END DO
         END DO

        IPERM=INDEX_RANK

#ifdef DEBUG_PRINT
        PRINT *,"*******  NO ENV NAG *****************"
        DO I=1,M2
          PRINT *,'IPERM(',I,')=',IPERM(I)
        END DO
#endif

#endif
      end subroutine MATH_QTLMAP_M01ZAF


!****************************************************************************************************************
!!  NAME
!!    IM01DFF
!!  DESCRIPTION
!!     ranks rows M1 to M2 of a matrix, using the data in columns N1 to N2 of those rows. The
!!    ordering is determined by first ranking the data in column N1, then ranking any tied rows according to the
!!    data in column N1 + 1, and so on up to column N2.


      subroutine MATH_QTLMAP_M01DFF(IM, LDM, M1, M2, N1, N2, ORDER, IRANK, IFAIL)
        INTEGER                     , INTENT(IN)    :: LDM,M1,M2,N1,N2
        INTEGER   ,DIMENSION(LDM,N2), INTENT(IN)    :: IM ! column of N1 to N2 of rows M1 to M2 of IM contains INTEGER to be ranked
        CHARACTER(LEN=1)            , INTENT(IN)    :: ORDER  ! 'A' ascending or 'D' descending order
        INTEGER   ,DIMENSION(M2), INTENT(OUT)       :: IRANK  ! elements M1 to M2 of IRANK contain the ranks of the corresponding rows of IM.
!  Note that the ranks are in the range M1 to M2: thus, if the ith row of IM is the first in the rank order,
!  IRANK(i) is set to M1.

        INTEGER                 , INTENT(OUT)       :: IFAIL
#ifdef HAVE_LIBNAG
          external :: M01DFF
#ifdef DEBUG_PRINT
            PRINT *," --> M01DFF ** nag **"
#endif
          call M01DFF(IM, LDM, M1, M2, N1, N2, ORDER, IRANK, IFAIL)
#ifdef DEBUG_PRINT
            PRINT *,"*********** ENV NAG ***********"
            DO I=M1,M2
              PRINT *,'IRANK(',I,')=',IRANK(I)
            END DO
#endif
#else
          INTEGER , DIMENSION(N1:N2,M2) :: INT_IRANK,IIRANK
          INTEGER                       :: I,J,K,IND,IND2,CURRENT_RANK,AFFECTED(M1:M2)
          LOGICAL                       :: FILTERIN(M1:M2),FILTEROUT(M1:M2)
#ifdef DEBUG_PRINT
            PRINT *," --> IM01DFF"
#endif

          IF ( (M2<1).or.(N2<1).or.(M1<1).or.(M1>M2).or.(N1<1).or.(N1>N2).or.(LDM<M2)) THEN
           IFAIL = 1
           RETURN
          END IF

          IF ( ORDER /= 'A' .and. ORDER /= 'a' .and. ORDER /= 'D' .and. ORDER /= 'd') THEN
           IFAIL = 2
           RETURN
          END IF

          ! Build rank vector for each column

          DO I=N1,N2
            CALL MATH_QTLMAP_M01DBF (IM(:,I), M1, M2, ORDER, INT_IRANK(I,:), IFAIL)
            IF ( IFAIL /= 0 ) THEN
               PRINT *,'IM01DFF-->IM01DBF IFAIL=',IFAIL
               RETURN
            END IF
#ifdef DEBUG_PRINT
            PRINT *,"IM01DBF WITH COLUMN:",I
            DO J=M1,M2
              PRINT *,'INT_IRANK(',J,')=',INT_IRANK(I,J)
            END DO
#endif
          END DO ! END I

         !Table of index rank
         DO I=N1,N2
           DO J=M1,M2
             IIRANK(I,INT_IRANK(I,J))=J
           END DO
         END DO



          IRANK=0
          FILTERIN=.true.
          FILTEROUT=.false.
          DO CURRENT_RANK=M1,M2

           DO J=N1,N2
              IF ( .not. is_conflict(LDM,M1,M2,N1,N2,J,FILTERIN,FILTEROUT,IM,INT_IRANK,IIRANK) .or. J==N2)  THEN
               !find the rows
               DO I=M1,M2
                IF ( FILTEROUT(I) ) THEN
                 EXIT
                END IF
               END DO
               !set the rank
               IRANK(I)=CURRENT_RANK
               AFFECTED(CURRENT_RANK)=I
               FILTERIN=.true.
               DO I=M1,CURRENT_RANK
                 FILTERIN(AFFECTED(I))=.false.
               END DO
               EXIT

              ELSE ! THE FILTER OUT contains several TRUE value, we have to search in the next row
                 !initializing the filterin
                 FILTERIN=FILTEROUT
              END IF
            END DO ! J
          END DO  !I

#ifdef DEBUG_PRINT
            PRINT *,"*********** NO ENV NAG ***********"
            DO I=M1,M2
              PRINT *,'IRANK(',I,')=',IRANK(I)
            END DO
#endif

#endif
      end subroutine MATH_QTLMAP_M01DFF

               ! Get the next bolean colum in the jth column from a boolean table (jth -1 column)
               ! 3 2 1  for the rank=1  ==> false  ....  ....
               ! 2 2 1                      true   true  true
               ! 2 2 2                      true   true  ...
               ! 2 3 1                      true   false  ...
      function is_conflict(LDM,M1,M2,N1,N2,J,tabtocomp,nexttocomp,IM,IIRANK,INDEXIRANK) result(RET)
            integer                   , INTENT(IN)      :: M1,M2,J,LDM,N1,N2
            logical,dimension(M1:M2)  , INTENT(IN)      :: tabtocomp
            logical,dimension(M1:M2)  , INTENT(OUT)     :: nexttocomp
            INTEGER   ,DIMENSION(LDM,N2), INTENT(IN)    :: IM
            INTEGER   ,DIMENSION(N1:N2,M2), INTENT(IN)     :: IIRANK,INDEXIRANK

            INTEGER :: I,IND,MINV,TABIND(count(tabtocomp)),SUBIIRANK(count(tabtocomp))
            LOGICAL :: RET
#ifdef DEBUG_PRINT
            PRINT *,'** CREATE SUBTAB ** SIZE:',SIZE(TABIND)
#endif
            IND=1
            DO I=M1,M2
              IF (tabtocomp(I)) THEN
#ifdef DEBUG_PRINT
                PRINT *,'ADD ROW:',I,' RANK INDEX:',IIRANK(J,I)
#endif
                TABIND(IND)=I
                SUBIIRANK(IND)=IIRANK(J,I)
                IND = IND+1
              END IF
            END DO

            !SEARCH THE MINIMUM RANK
            minv=minval(SUBIIRANK)
#ifdef DEBUG_PRINT
            PRINT *,'RANK VALMIN:',minv,' ROWS:',INDEXIRANK(J,minv),' corresponding to the value:',IM(INDEXIRANK(J,minv),J)
#endif
            !CREATE NEW TAB FILTER BOOLEAN
            nexttocomp=.false.
            DO I=1,SIZE(TABIND)
              IND=TABIND(I)
              IF ( IM(IND,J) == IM(INDEXIRANK(J,minv),J) ) THEN
#ifdef DEBUG_PRINT
                 PRINT *,'-->KEEP ROW:',IND
#endif
                 nexttocomp(IND)=.true.
              END IF
            END DO

            !IF ONE VALUE=>NO CONFLICT,WE FIND THE ROW
            IF ( count(nexttocomp) == 1 ) THEN
               RET=.false.
            ELSE
               RET=.true.
            END IF

            RETURN

      end function is_conflict

!****************************************************************************************************************

!!****f* m_qtlmap_math/MATH_QTLMAP_F03ABF
!!  NAME
!!    MATH_QTLMAP_F03ABF
!!  DESCRIPTION
!!    Calculates the determinant of a real symmetric positive-definite matrix using a Cholesky factorization.
!!  NOTES
!!
!!***
      subroutine MATH_QTLMAP_F03ABF(A, IA, N, DET, IFAIL)
          integer      ,intent(in)                   :: IA,N
          double precision,dimension(IA,N), intent(in)  :: A
          double precision,intent(out)                  :: DET
          integer      ,intent(out)                  :: IFAIL

#ifdef HAVE_LIBNAG
          double precision,dimension(N)     :: WKSPCE
#else
          double precision,dimension(IA,N)  :: APRIM
          double precision,dimension(N)     :: eigenv
          double precision                  :: trace, sume
          double precision,parameter        :: TOL = 0.001d0

          IFAIL = 0
#endif
          DET=0.d0
#ifdef HAVE_LIBNAG
          call F03ABF(A, IA, N, DET, WKSPCE, IFAIL)
#else

          APRIM=A
          call kaiser(APRIM, IA, N, eigenv, trace, sume, IFAIL)

          ! Probleme de precision....
          if (trace == 0.d0) then
             IFAIL=3
             DET=0.d0
             RETURN
          end if
          IF ( abs((trace-sume)/trace)>TOL ) THEN
             DET = product(eigenv)
             IFAIL= 3
          END IF

          IF ( IFAIL == 0 ) THEN
            DET = product(eigenv)
          ELSE
            DET = 0.d0
          END IF

#endif
#ifdef DEBUG
          call log_mess('MATH_QTLMAP_F03ABF (OUT) ==== DET:'//str(DET),DEBUG_DEF)
#endif
      end subroutine MATH_QTLMAP_F03ABF

!****************************************************************************************************************
!!****f* m_qtlmap_math/MATH_QTLMAP_F01ADF
!!  NAME
!!    MATH_QTLMAP_F01ADF
!!  DESCRIPTION
!!   Calculates the approximate inverse of a real symmetric positive-definite matrix, using a Cholesky
!!   factorization.
!!  NOTES
!!   USe DPOTRF      : LAPACK
!!
!!***
       subroutine MATH_QTLMAP_F01ADF (N, A, IA, IFAIL)
          integer         ,intent(in)                     :: IA,N
          double precision,dimension(IA,N) ,intent(inout) :: A
          integer         ,intent(out)                    :: IFAIL
#ifdef HAVE_LIBNAG

#else
          double precision,dimension(IA,N)  :: APRIM
          double precision,dimension(2)     :: DET2
          integer                           :: INFO,I,J
         ! external                          :: DPOFA,DPODI ! SLATEC
          external                          :: DPOTRF,DPODI ! LAPACK
          character(len=1)                  :: UPLO
          IFAIL=0
#endif

#ifdef HAVE_LIBNAG
         call F01ADF(N, A, IA, IFAIL)
#else
          UPLO = 'U'
          APRIM(1:N,1:N)=A(1:N,1:N)
          ! CALL DPOFA(APRIM, IA, N, INFO)   ! ** SLATEC is DEPRECATED **
          CALL DPOTRF( UPLO, N, APRIM, IA, INFO )

          IF ( INFO /= 0 ) THEN
              IFAIL=1
              RETURN
          END IF

          !CALL DPODI(APRIM,IA,N,DET2,01) ! ** SLATEC is DEPRECATED **

          CALL DPOTRI( UPLO, N, APRIM, IA, INFO )

          IF ( INFO /= 0 ) THEN
              IFAIL=1
              RETURN
          END IF

          DO I=2,IA
            DO J=1,(I-1)
                A(I,J)=APRIM(J,(I-1))
            END DO
          END DO

!          print *,"**MATH_QTLMAP_F01ADF **"
!          DO I=1,I
!            PRINT *,A(I,:)
!          END DO
#endif
       end subroutine MATH_QTLMAP_F01ADF

!****************************************************************************************************************
!!****f* m_qtlmap_math/MATH_QTLMAP_INVDETMATSYM
!!  NAME
!!    MATH_QTLMAP_F01ADF
!!  DESCRIPTION
!!   Calculates the inverse and the derterminant of a real symmetric positive-definite matrix
!!  NOTES
!!
!!***
       subroutine MATH_QTLMAP_INVDETMATSYM (N, A, IA,DET, IFAIL)
          integer         ,intent(in)                     :: IA,N
          double precision,dimension(IA,N) ,intent(inout) :: A
          double precision                   ,intent(out) :: DET
          integer         ,intent(out)                    :: IFAIL

          double precision,dimension(IA,N)  :: APRIM
          double precision,dimension(2)     :: DET2
          integer                           :: INFO,I,J
          external                          :: DPOTRF,DPODI ! LAPACK
          character(len=1)                  :: UPLO
          IFAIL=0

          UPLO = 'U'
          APRIM(1:N,1:N)=A(1:N,1:N)
          ! CALL DPOFA(APRIM, IA, N, INFO)   ! ** SLATEC is DEPRECATED **
          CALL DPOTRF( UPLO, N, APRIM, IA, INFO )

          IF ( INFO /= 0 ) THEN
              IFAIL=1
              RETURN
          END IF
          DET=1
          do i=1,N
            DET=DET*APRIM(I,I)
          end do
          DET=DET**2

          CALL DPOTRI( UPLO, N, APRIM, IA, INFO )

          IF ( INFO /= 0 ) THEN
              IFAIL=1
              RETURN
          END IF

          DO I=2,IA
            DO J=1,(I-1)
                A(I,J)=APRIM(J,(I-1))
            END DO
          END DO

       end subroutine MATH_QTLMAP_INVDETMATSYM


!!****f* m_qtlmap_math/MATH_QTLMAP_INVDETMAT
!!  NAME
!!    MATH_QTLMAP_INVDETMAT
!!  DESCRIPTION
!!   Calculates the inverse and the derterminant of a real matrix
!!  NOTES
!!   Use DPOFA,DPODI : SLATEC
!!***
!       subroutine MATH_QTLMAP_INVDETMAT (N, A, DET, IFAIL)
!          integer         ,intent(in)                     :: N
!          double precision,dimension(N,N) ,intent(inout)  :: A
!          double precision                   ,intent(out) :: DET
!          integer         ,intent(out)                    :: IFAIL
!
!          double precision,dimension(2)     :: DET2
!          integer                           :: INFO,I,J
!          integer         ,dimension(N)     :: KPVT
!          integer         ,dimension(3)     :: INERT
!          external                          :: DSIFA,DSIDI ! SLATEC
!          double precision,dimension(N)     :: WORK
!          integer                           :: JOB
!          IFAIL=0
!
!          CALL DSIFA(A, N, N, KPVT,INFO)
!          IF ( INFO /= 0 ) THEN
!              IFAIL=1
!              write (0,*)"ERROR:",INFO
!              CALL log_mess('DSIFA --> Can not compute inverse of matrix',ERROR_DEF)
!              RETURN
!          END IF
!          JOB=110
!          CALL DSIDI(A,N,N,KPVT,DET2,INERT,WORK,JOB)
!          DET = DET2(1) * 10.d0**DET2(2)
!
!       end subroutine MATH_QTLMAP_INVDETMAT


!!****f* m_qtlmap_math/MATH_QTLMAP_S15AEF
!!  NAME
!!    MATH_QTLMAP_S15AEF
!!  DESCRIPTION
!!   Returns the value of the error function erf x, via the routine name.
!!  NOTES
!!   Use DERF from http://www.kurims.kyoto-u.ac.jp/~ooura/gamerf.html
!!***
    function MATH_QTLMAP_S15AEF(X,IFAIL) result(y)
        double precision,intent(in)  :: X
        integer, intent(inout)       :: IFAIL
        double precision             :: y

#ifdef HAVE_LIBNAG
        double precision,external    :: S15AEF
         y = S15AEF(X, IFAIL)
#else
         y = derf(X)
#endif
    end function MATH_QTLMAP_S15AEF
!****************************************************************************************************************
!!****f* m_qtlmap_math/MATH_QTLMAP_S15ADF
!!  NAME
!!    MATH_QTLMAP_S15ADF
!!  DESCRIPTION
!!   Returns the value of the complementary error function, erfc x, via the routine name.
!!  NOTES
!!   Use DERFC from http://www.kurims.kyoto-u.ac.jp/~ooura/gamerf.html
!!***
    function MATH_QTLMAP_S15ADF(X,IFAIL) result(y)
        double precision,intent(in)  :: X
        integer, intent(inout)       :: IFAIL
        double precision             :: y
#ifdef HAVE_LIBNAG
        double precision,external    :: S15ADF
         y = S15ADF(X, IFAIL)
#else
         y = derfc(X)
#endif
    end function MATH_QTLMAP_S15ADF

!****************************************************************************************************************

       function MATH_QTLMAP_GAMMA_F(N) result(R)
            integer,intent(in)   :: N
            integer :: P,R

            P = N-1
            R = P
            do while ( P > 2 )
              P = P - 1
              R = R*P
            end do

       end function MATH_QTLMAP_GAMMA_F

! ren-raw chen, rutgers business school
! normal inverse
! translate from
!http://home.online.no/~pjacklam/notes/invnorm
! a routine written by john herrero
      function dinvnorm(p) result (X)
      double precision,intent(in) :: p
      double precision :: p_low,p_high
      double precision :: a1,a2,a3,a4,a5,a6
      double precision :: b1,b2,b3,b4,b5
      double precision :: c1,c2,c3,c4,c5,c6
      double precision :: d1,d2,d3,d4
      double precision :: z,q,r,X

      a1=-39.6968302866538
      a2=220.946098424521
      a3=-275.928510446969
      a4=138.357751867269
      a5=-30.6647980661472
      a6=2.50662827745924
      b1=-54.4760987982241
      b2=161.585836858041
      b3=-155.698979859887
      b4=66.8013118877197
      b5=-13.2806815528857
      c1=-0.00778489400243029
      c2=-0.322396458041136
      c3=-2.40075827716184
      c4=-2.54973253934373
      c5=4.37466414146497
      c6=2.93816398269878
      d1=0.00778469570904146
      d2=0.32246712907004
      d3=2.445134137143
      d4=3.75440866190742
      p_low=0.02425
      p_high=1-p_low
      if(p.lt.p_low) goto 201
      if(p.ge.p_low) goto 301
201   q=dsqrt(-2*dlog(p))
      z=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/ &
     ((((d1*q+d2)*q+d3)*q+d4)*q+1)
      goto 204
301   if((p.ge.p_low).and.(p.le.p_high)) goto 202
      if(p.gt.p_high) goto 302
202   q=p-0.5
      r=q*q
      z=(((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q/ &
     (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1)
      goto 204
302   if((p.gt.p_high).and.(p.lt.1)) goto 203
203   q=dsqrt(-2*dlog(1-p))
      z=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/  &
     ((((d1*q+d2)*q+d3)*q+d4)*q+1)
204   X=z
      return
      end function dinvnorm


! error function in double precision
! http://www.kurims.kyoto-u.ac.jp/~ooura/gamerf.html
      function derf(x) result(P)
       double precision  ,intent(in)      :: x
       double precision , dimension(0:64) :: a, b
       double precision                   :: w,t,y,P
       integer                            :: i,k

      data (a(i), i = 0, 12) /                                 &
         0.00000000005958930743d0, -0.00000000113739022964d0,  &
         0.00000001466005199839d0, -0.00000016350354461960d0,  &
         0.00000164610044809620d0, -0.00001492559551950604d0,  &
         0.00012055331122299265d0, -0.00085483269811296660d0,  &
         0.00522397762482322257d0, -0.02686617064507733420d0,  &
         0.11283791670954881569d0, -0.37612638903183748117d0,  &
         1.12837916709551257377d0 /
      data (a(i), i = 13, 25) /                                &
         0.00000000002372510631d0, -0.00000000045493253732d0,  &
         0.00000000590362766598d0, -0.00000006642090827576d0,  &
         0.00000067595634268133d0, -0.00000621188515924000d0,  &
         0.00005103883009709690d0, -0.00037015410692956173d0,  &
         0.00233307631218880978d0, -0.01254988477182192210d0,  &
         0.05657061146827041994d0, -0.21379664776456006580d0,  &
         0.84270079294971486929d0 /
      data (a(i), i = 26, 38) /                                &
         0.00000000000949905026d0, -0.00000000018310229805d0,  &
         0.00000000239463074000d0, -0.00000002721444369609d0,  &
         0.00000028045522331686d0, -0.00000261830022482897d0,  &
         0.00002195455056768781d0, -0.00016358986921372656d0,  &
         0.00107052153564110318d0, -0.00608284718113590151d0,  &
         0.02986978465246258244d0, -0.13055593046562267625d0,  &
         0.67493323603965504676d0 /
      data (a(i), i = 39, 51) /                                &
         0.00000000000382722073d0, -0.00000000007421598602d0,  &
         0.00000000097930574080d0, -0.00000001126008898854d0,  &
         0.00000011775134830784d0, -0.00000111992758382650d0,  &
         0.00000962023443095201d0, -0.00007404402135070773d0,  &
         0.00050689993654144881d0, -0.00307553051439272889d0,  &
         0.01668977892553165586d0, -0.08548534594781312114d0,  &
         0.56909076642393639985d0 /
      data (a(i), i = 52, 64) /                                &
         0.00000000000155296588d0, -0.00000000003032205868d0,  &
         0.00000000040424830707d0, -0.00000000471135111493d0,  &
         0.00000005011915876293d0, -0.00000048722516178974d0,  &
         0.00000430683284629395d0, -0.00003445026145385764d0,  &
         0.00024879276133931664d0, -0.00162940941748079288d0,  &
         0.00988786373932350462d0, -0.05962426839442303805d0,  &
         0.49766113250947636708d0 /
      data (b(i), i = 0, 12) /                                 &
         -0.00000000029734388465d0, 0.00000000269776334046d0,  &
         -0.00000000640788827665d0, -0.00000001667820132100d0, &
         -0.00000021854388148686d0, 0.00000266246030457984d0,  &
         0.00001612722157047886d0, -0.00025616361025506629d0,  &
         0.00015380842432375365d0, 0.00815533022524927908d0,  &
         -0.01402283663896319337d0, -0.19746892495383021487d0,  &
         0.71511720328842845913d0 /
      data (b(i), i = 13, 25) /                                  &
         -0.00000000001951073787d0, -0.00000000032302692214d0,   &
         0.00000000522461866919d0, 0.00000000342940918551d0,     &
         -0.00000035772874310272d0, 0.00000019999935792654d0,    &
         0.00002687044575042908d0, -0.00011843240273775776d0,    &
         -0.00080991728956032271d0, 0.00661062970502241174d0,    &
         0.00909530922354827295d0, -0.20160072778491013140d0,    &
         0.51169696718727644908d0 /
      data (b(i), i = 26, 38) /                                  &
         0.00000000003147682272d0, -0.00000000048465972408d0,    &
         0.00000000063675740242d0, 0.00000003377623323271d0,     &
         -0.00000015451139637086d0, -0.00000203340624738438d0,   &
         0.00001947204525295057d0, 0.00002854147231653228d0,     &
         -0.00101565063152200272d0, 0.00271187003520095655d0,    &
         0.02328095035422810727d0, -0.16725021123116877197d0,    &
         0.32490054966649436974d0 /
      data (b(i), i = 39, 51) /                                  &
         0.00000000002319363370d0, -0.00000000006303206648d0,    &
         -0.00000000264888267434d0, 0.00000002050708040581d0,    &
         0.00000011371857327578d0, -0.00000211211337219663d0,    &
         0.00000368797328322935d0, 0.00009823686253424796d0,     &
         -0.00065860243990455368d0, -0.00075285814895230877d0,   &
         0.02585434424202960464d0, -0.11637092784486193258d0,    &
         0.18267336775296612024d0 /
      data (b(i), i = 52, 64) /                                  &
         -0.00000000000367789363d0, 0.00000000020876046746d0,    &
         -0.00000000193319027226d0, -0.00000000435953392472d0,   &
         0.00000018006992266137d0, -0.00000078441223763969d0,    &
         -0.00000675407647949153d0, 0.00008428418334440096d0,    &
         -0.00017604388937031815d0, -0.00239729611435071610d0,   &
         0.02064129023876022970d0, -0.06905562880005864105d0,    &
         0.09084526782065478489d0 /
      w = abs(x)
      if (w .lt. 2.2d0) then
          t = w * w
          k = int(t)
          t = t - k
          k = k * 13
          y = ((((((((((((a(k) * t + a(k + 1)) * t +              &
             a(k + 2)) * t + a(k + 3)) * t + a(k + 4)) * t +      &
             a(k + 5)) * t + a(k + 6)) * t + a(k + 7)) * t +       &
             a(k + 8)) * t + a(k + 9)) * t + a(k + 10)) * t +     &
             a(k + 11)) * t + a(k + 12)) * w
      else if (w .lt. 6.9d0) then
          k = int(w)
          t = w - k
          k = 13 * (k - 2)
          y = (((((((((((b(k) * t + b(k + 1)) * t + &
             b(k + 2)) * t + b(k + 3)) * t + b(k + 4)) * t + &
             b(k + 5)) * t + b(k + 6)) * t + b(k + 7)) * t + &
             b(k + 8)) * t + b(k + 9)) * t + b(k + 10)) * t + &
             b(k + 11)) * t + b(k + 12)
          y = y * y
          y = y * y
          y = y * y
          y = 1 - y * y
      else
          y = 1
      end if
      if (x .lt. 0) y = -y
      P = y
      end function derf
!

      ! complementary error function in double precision
!! http://www.kurims.kyoto-u.ac.jp/~ooura/gamerf.html
      function derfc(x) result(y)
      double precision , intent(in) :: x
      double precision              :: y,t,u

      double precision,parameter :: pa = 3.97886080735226000d+00
         double precision,parameter ::  p0 = 2.75374741597376782d-01
         double precision,parameter ::  p1 = 4.90165080585318424d-01
         double precision,parameter ::  p2 = 7.74368199119538609d-01
         double precision,parameter ::  p3 = 1.07925515155856677d+00
         double precision,parameter ::  p4 = 1.31314653831023098d+00
         double precision,parameter ::  p5 = 1.37040217682338167d+00
         double precision,parameter ::  p6 = 1.18902982909273333d+00
         double precision,parameter ::  p7 = 8.05276408752910567d-01
         double precision,parameter ::  p8 = 3.57524274449531043d-01
         double precision,parameter ::  p9 = 1.66207924969367356d-02
         double precision,parameter ::  p10 = -1.19463959964325415d-01
         double precision,parameter ::  p11 = -8.38864557023001992d-02

         double precision,parameter ::  p12 = 2.49367200053503304d-03
         double precision,parameter ::  p13 = 3.90976845588484035d-02
         double precision,parameter ::  p14 = 1.61315329733252248d-02
         double precision,parameter ::  p15 = -1.33823644533460069d-02
         double precision,parameter ::  p16 = -1.27223813782122755d-02
         double precision,parameter ::  p17 = 3.83335126264887303d-03
         double precision,parameter ::  p18 = 7.73672528313526668d-03
         double precision,parameter ::  p19 = -8.70779635317295828d-04
         double precision,parameter ::  p20 = -3.96385097360513500d-03
         double precision,parameter ::  p21 = 1.19314022838340944d-04
         double precision,parameter ::  p22 = 1.27109764952614092d-03

      t = pa / (pa + abs(x))
      u = t - 0.5d0
      y = (((((((((p22 * u + p21) * u + p20) * u +         &
         p19) * u + p18) * u + p17) * u + p16) * u +       &
         p15) * u + p14) * u + p13) * u + p12
      y = ((((((((((((y * u + p11) * u + p10) * u +        &
         p9) * u + p8) * u + p7) * u + p6) * u + p5) * u + &
         p4) * u + p3) * u + p2) * u + p1) * u + p0) * t * &
         exp(-x * x)
      if (x .lt. 0) y = 2 - y

      end function derfc
!

! Gamma function in double precision
! ! http://www.kurims.kyoto-u.ac.jp/~ooura/gamerf.html
      function dgamma(x) result(RES)
      double precision , intent(in) :: x
      double precision              :: w,y,RES
      integer                       :: k,n

         double precision,parameter :: p0 = 0.999999999999999990d+00
         double precision,parameter :: p1 = -0.422784335098466784d+00
         double precision,parameter :: p2 = -0.233093736421782878d+00
         double precision,parameter :: p3 = 0.191091101387638410d+00
         double precision,parameter :: p4 = -0.024552490005641278d+00
         double precision,parameter :: p5 = -0.017645244547851414d+00
         double precision,parameter :: p6 = 0.008023273027855346d+00

         double precision,parameter :: p7 = -0.000804329819255744d+00
         double precision,parameter :: p8 = -0.000360837876648255d+00
         double precision,parameter :: p9 = 0.000145596568617526d+00
         double precision,parameter :: p10 = -0.000017545539395205d+00
         double precision,parameter :: p11 = -0.000002591225267689d+00
         double precision,parameter :: p12 = 0.000001337767384067d+00
         double precision,parameter :: p13 = -0.000000199542863674d+00
      n = nint(x - 2)
      w = x - (n + 2)
      y = ((((((((((((p13 * w + p12) * w + p11) * w + p10) *  &
         w + p9) * w + p8) * w + p7) * w + p6) * w + p5) *    &
         w + p4) * w + p3) * w + p2) * w + p1) * w + p0
      if (n .gt. 0) then
          w = x - 1
          do k = 2, n
              w = w * (x - k)
          end do
      else
          w = 1
          do k = 0, -n - 1
              y = y * (x + k)
          end do
      end if
      RES = w / y
      end function dgamma


      SUBROUTINE kaiser(a, nrows, n, eigenv, trace, sume, ier)

!  EIGENVALUES AND VECTORS OF A SYMMETRIC +VE DEFINITE MATRIX,
!  USING KAISER'S METHOD.
!  REFERENCE: KAISER,H.F. 'THE JK METHOD: A PROCEDURE FOR FINDING THE
!  EIGENVALUES OF A REAL SYMMETRIC MATRIX', COMPUT.J., VOL.15, 271-273, 1972.

!  ARGUMENTS:-
!  A       = INPUT, AN ARRAY CONTAINING THE MATRIX
!            OUTPUT, THE COLUMNS OF A CONTAIN THE NORMALIZED EIGENVECTORS
!            OF A.   N.B. A IS OVERWRITTEN !
!  NROWS   = INPUT, THE FIRST DIMENSION OF A IN THE CALLING PROGRAM.
!  N       = INPUT, THE ORDER OF A, I.E. NO. OF COLUMNS.
!            N MUST BE <= NROWS.
!  EIGENV()= OUTPUT, A VECTOR CONTAINING THE ORDERED EIGENVALUES.
!  TRACE   = OUTPUT, THE TRACE OF THE INPUT MATRIX.
!  SUME    = OUTPUT, THE SUM OF THE EIGENVALUES COMPUTED.
!            N.B. ANY SYMMETRIC MATRIX MAY BE INPUT, BUT IF IT IS NOT +VE
!            DEFINITE, THE ABSOLUTE VALUES OF THE EIGENVALUES WILL BE FOUND.
!            IF TRACE = SUME, THEN ALL OF THE EIGENVALUES ARE POSITIVE
!            OR ZERO.   IF SUME > TRACE, THE DIFFERENCE IS TWICE THE SUM OF
!            THE EIGENVALUES WHICH HAVE BEEN GIVEN THE WRONG SIGNS !
!  IER     = OUTPUT, ERROR INDICATOR
!             = 0 NO ERROR
!             = 1 N > NROWS OR N < 1
!             = 2 FAILED TO CONVERGE IN 10 ITERATIONS

!  LATEST REVISION - 6 September 1990
!  Fortran 90 version - 20 November 1998

!*************************************************************************


REAL (kind=dp), INTENT(IN OUT) :: a(:,:)
INTEGER, INTENT(IN)       :: nrows
INTEGER, INTENT(IN)       :: n
REAL (kind=dp), INTENT(OUT)    :: eigenv(:)
REAL (kind=dp), INTENT(OUT)    :: trace
REAL (kind=dp), INTENT(OUT)    :: sume
INTEGER, INTENT(OUT)      :: ier

! Local variables

REAL (kind=dp), PARAMETER :: small = 1.0e-12_dp, zero = 0.0_dp, half = 0.5_dp,  &
                        one = 1.0_dp
INTEGER   :: i, iter, j, k, ncount, nn
REAL (kind=dp) :: absp, absq, COS, ctn, eps, halfp, p, q, SIN, ss, TAN, temp, xj, xk

!   CALCULATE CONVERGENCE TOLERANCE, EPS.
!   CALCULATE TRACE.   INITIAL SETTINGS.

ier = 1
IF(n < 1 .OR. n > nrows) RETURN
ier = 0
iter = 0
trace = zero
ss = zero
DO j = 1,n
  trace = trace + a(j,j)
  DO i = 1,n
    ss = ss + a(i,j)**2
  END DO
END DO
sume = zero
eps = small*ss/n
nn = n*(n-1)/2
ncount = nn

!   ORTHOGONALIZE PAIRS OF COLUMNS J & K, K > J.

20 DO j = 1,n-1
  DO k = j+1,n

!   CALCULATE PLANAR ROTATION REQUIRED

    halfp = zero
    q = zero
    DO i = 1,n
      xj = a(i,j)
      xk = a(i,k)
      halfp = halfp + xj*xk
      q = q + (xj+xk) * (xj-xk)
    END DO
    p = halfp + halfp
    absp = ABS(p)

!   If P is very small, the vectors are almost orthogonal.
!   Skip the rotation if Q >= 0 (correct ordering).

    IF (absp < eps .AND. q >= zero) THEN
      ncount = ncount - 1
      IF (ncount <= 0) GO TO 160
      CYCLE
    END IF

!   Rotation needed.

    absq = ABS(q)
    IF(absp <= absq) THEN
      TAN = absp/absq
      COS = one/SQRT(one + TAN*TAN)
      SIN = TAN*COS
    ELSE
      ctn = absq/absp
      SIN = one/SQRT(one + ctn*ctn)
      COS = ctn*SIN
    END IF
    COS = SQRT((one + COS)*half)
    SIN = SIN/(COS + COS)
    IF(q < zero) THEN
      temp = COS
      COS = SIN
      SIN = temp
    END IF
    IF(p < zero) SIN = -SIN

!   PERFORM ROTATION

    DO i = 1,n
      temp = a(i,j)
      a(i,j) = temp*COS + a(i,k)*SIN
      a(i,k) = -temp*SIN + a(i,k)*COS
    END DO
  END DO
END DO
ncount = nn
iter = iter + 1
IF(iter < 10) GO TO 20
ier = 2

!   CONVERGED, OR GAVE UP AFTER 10 ITERATIONS

160 DO j = 1,n
  temp = SUM( a(1:n,j)**2 )
  eigenv(j) = SQRT(temp)
  sume = sume + eigenv(j)
END DO

!   SCALE COLUMNS TO HAVE UNIT LENGTH

DO j = 1,n
  IF (eigenv(j) > zero) THEN
    temp = one/eigenv(j)
  ELSE
    temp = zero
  END IF
  a(1:n,j) = a(1:n,j)*temp
END DO

RETURN
END SUBROUTINE kaiser

!Student's t-distribution. (17/1968)
!START OF AS 3
      FUNCTION PROBST(T, IDF, IFAULT) result(prob)
!
!        ALGORITHM AS 3  APPL. STATIST. (1968) VOL.17, P.189
!
!        STUDENT T PROBABILITY (LOWER TAIL)
!
      DOUBLE PRECISION , INTENT(IN) :: T
      INTEGER          , INTENT(IN) :: IDF
      INTEGER          ,INTENT(OUT) :: IFAULT
      DOUBLE PRECISION :: A, B, C, F, G1, S, FK, ZERO, ONE, TWO, HALF, ZSQRT, ZATAN
      DOUBLE PRECISION :: prob
      INTEGER :: IM2,IOE,K,KS
!
!        G1 IS RECIPROCAL OF PI
!
      DATA ZERO, ONE, TWO, HALF, G1/0.0, 1.0, 2.0,  0.5, 0.3183098862/
!
      ZSQRT(A) = SQRT(A)
      ZATAN(A) = ATAN(A)
!
      IFAULT = 1
      PROB = ZERO
      IF (IDF .LT. 1) RETURN
      IFAULT = 0
      F = IDF
      A = T / ZSQRT(F)
      B = F / (F + T ** 2)
      IM2 = IDF - 2
      IOE = MOD(IDF, 2)
      S = ONE
      C = ONE
      F = ONE
      KS = 2 + IOE
      FK = KS
      IF (IM2 .LT. 2) GOTO 20
      DO 10 K = KS, IM2, 2
      C = C * B * (FK - ONE) / FK
      S = S + C
      IF (S .EQ. F) GOTO 20
      F = S
      FK = FK + TWO
   10 CONTINUE
   20 IF (IOE .EQ. 1) GOTO 30
      PROB = HALF + HALF * A * ZSQRT(B) * S
      RETURN
   30 IF (IDF .EQ. 1) S = ZERO
      PROB = HALF + (A * B * S + ZATAN(A)) * G1
      RETURN
      END FUNCTION PROBST
!END OF AS 3

      FUNCTION STUDNT (T, DOFF, IFAULT) result(ut)
!
!     ALGORITHM AS 27  APPL. STATIST. VOL.19, NO.1
!
!     Calculate the upper tail area under Student's t-distribution
!
!     Translated from Algol by Alan Miller
!
      INTEGER IFAULT
      DOUBLE PRECISION T, DOFF,UT
!
!     Local variables
!
      DOUBLE PRECISION V, X, TT, TWO, FOUR, ONE, ZERO, HALF
      DOUBLE PRECISION A1, A2, A3, A4, A5, B1, B2, &
          C1, C2, C3, C4, C5, D1, D2,  &
          E1, E2, E3, E4, E5, F1, F2,  &
          G1, G2, G3, G4, G5, H1, H2,  &
          I1, I2, I3, I4, I5, J1, J2
      LOGICAL POS
      DATA TWO /2.0/, FOUR /4.0/, ONE /1.0/, ZERO /0.0/, HALF /0.5/
      DATA A1, A2, A3, A4, A5 /0.09979441, -0.581821, 1.390993,      &
          -1.222452, 2.151185/, B1, B2 /5.537409, 11.42343/
      DATA C1, C2, C3, C4, C5 /0.04431742, -0.2206018, -0.03317253,  &
          5.679969, -12.96519/, D1, D2 /5.166733, 13.49862/
      DATA E1, E2, E3, E4, E5 /0.009694901, -0.1408854, 1.88993,     &
          -12.75532, 25.77532/, F1, F2 /4.233736, 14.3963/
      DATA G1, G2, G3, G4, G5 /-9.187228E-5, 0.03789901, -1.280346,  &
          9.249528, -19.08115/, H1, H2 /2.777816, 16.46132/
      DATA I1, I2, I3, I4, I5 /5.79602E-4, -0.02763334, 0.4517029,   &
          -2.657697, 5.127212/, J1, J2 /0.5657187, 21.83269/
!
!     Check that number of degrees of freedom > 4.
!
      IF (DOFF .LT. TWO) THEN
	IFAULT = 1
	UT = - ONE
	RETURN
      END IF
!
      IF (DOFF .LE. FOUR) THEN
	IFAULT = DOFF
      ELSE
	IFAULT = 0
      END IF
!
!     Evaluate series.
!
      V = ONE / DOFF
      POS = (T .GE. ZERO)
      TT = ABS(T)
      X = HALF * (ONE +                                          &
         TT * (((A1 + V * (A2 + V * (A3 + V * (A4 + V * A5)))) / &
             (ONE - V * (B1 - V * B2))) +                        &
         TT * (((C1 + V * (C2 + V * (C3 + V * (C4 + V * C5)))) / &
             (ONE - V * (D1 - V * D2))) +                        &
         TT * (((E1 + V * (E2 + V * (E3 + V * (E4 + V * E5)))) / &
             (ONE - V * (F1 - V * F2))) +                        &
         TT * (((G1 + V * (G2 + V * (G3 + V * (G4 + V * G5)))) / &
             (ONE - V * (H1 - V * H2))) +                        &
         TT * ((I1 + V * (I2 + V * (I3 + V * (I4 + V * I5)))) /  &
             (ONE - V * (J1 - V * J2))) ))))) ** (-8)
      IF (POS) THEN
	UT = X
      ELSE
	UT = ONE - X
      END IF
!
      RETURN
      END FUNCTION STUDNT

!Complete beta function
  function betaC(x,y) result(RES)
    double precision , intent(in)  :: x,y
    double precision   :: RES

    RES = dgamma(x)*dgamma(y)/dgamma(x+Y)

  end function betaC


!Incomplete gamma function
      function gammds (y,p,ifault) result(RES)
!
!        Algorithm AS 147  Appl. Statist. (1980) Vol. 29, No. 1
!
!        Computes the incomplete gamma integral for positive
!        parameters y,p using an infinite series
!
!        Auxiliary function required: ALNGAM = CACM algorithm 291
!
!	     AS239 should be considered as an alternative to AS147
!
      double precision ,intent(in):: y,p
      integer, intent(out)  :: ifault
      double precision :: RES,e,zero,one,uflo,arg,c,a,f
      data e/1.0d-9/, zero/0.0d0/, one/1.0d0/, uflo/1.0d-37/
!
!        Checks admissibility of arguments and value of f
!
      ifault = 1
      RES = zero
      if(y.le.zero .or. p.le.zero) return
      ifault = 2
!
!        alngam is natural log of gamma function
!
      arg = p*log(y)-log(dgamma(p+one))-y
      if(arg.lt.log(uflo)) return
      f = exp(arg)
      if(f.eq.zero) return
      ifault = 0
!
!          Series begins
!
      c = one
      RES = one
      a = p
    1 a = a+one
      c = c*y/a
      RES = RES+c
      if (c/RES.gt.e) goto 1
      RES = RES*f
      return
      end function gammds




FUNCTION betain(x, p, q, beta) RESULT(fn_val)

!   Algorithm AS 63  Appl. Statist. (1973), vol.22, no.3

!   Computes incomplete beta function ratio for arguments
!   x between zero and one, p and q positive.
!   Log of complete beta function, beta, is assumed to be known

! ELF90-compatible version by Alan Miller
! N.B. Argument IFAULT has been removed

! Latest revision - 5 July 2003

USE LogBeta
IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(15, 100)
REAL (dp), INTENT(IN) :: x, p, q, beta
REAL (dp)             :: fn_val

! Local variables
LOGICAL               :: indx
INTEGER               :: ns
REAL (dp)             :: psq, cx, xx, pp, qq, term, ai, rx, temp

!   Define accuracy and initialise

REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp, acu = 1.0E-14_dp

fn_val = x

!   Test for admissibility of arguments

IF(p <= zero .OR. q <= zero) THEN
  WRITE(*, *) 'AS63: Either p or q <= 0'
  RETURN
END IF
IF(x < zero .OR. x > one) THEN
  WRITE(*, *) 'AS63: Argument x outside range (0, 1)'
  RETURN
END IF
IF(x == zero .OR. x == one) RETURN

!   Change tail if necessary and determine s

psq = p + q
cx = one - x
IF(p < psq*x) THEN
  xx = cx
  cx = x
  pp = q
  qq = p
  indx = .true.
ELSE
  xx = x
  pp = p
  qq = q
  indx = .false.
END IF
term = one
ai = one
fn_val = one
ns = qq + cx*psq

!     Use Soper's reduction formulae.

rx = xx/cx
3 temp = qq - ai
IF(ns == 0) rx = xx
4 term = term*temp*rx/(pp+ai)
fn_val = fn_val + term
temp = ABS(term)
IF(temp <= acu .AND. temp <= acu*fn_val) GO TO 5
ai = ai + one
ns = ns - 1
IF(ns >= 0) GO TO 3
temp = psq
psq = psq+one
GO TO 4

!     Calculate result

5 fn_val = fn_val*EXP(pp*LOG(xx) + (qq-one)*LOG(cx) - beta)/pp
IF(indx) fn_val = one - fn_val

RETURN
END FUNCTION betain

      FUNCTION DBETAI (X, PIN, QIN) RESULT (RES)
         use LogBeta
         REAL (dp) :: X, PIN, QIN,beta,RES



         beta = betaln(PIN, QIN)
         RES = betain(X, PIN, QIN, beta)
         !print *,"DBETAI:::::",RES


      END FUNCTION DBETAI

!!SLATEC
!!------

!     FUNCTION D1MACH (I) RESULT(RES)
!         INTEGER , INTENT(IN) :: I
!         DOUBLE PRECISION :: RES
!         IF ( I == 1 ) RES = tiny ( RES )
!         IF ( I == 2 ) RES = huge ( RES )
!         IF ( I == 3 ) RES = epsilon ( RES ) / radix ( RES )
!         IF ( I == 4 ) RES = epsilon ( RES )
!         IF ( I == 5 ) RES = log ( dble(radix ( RES )) )
!     END FUNCTION D1MACH


!*DECK DBETAI
 !     FUNCTION DBETAI (X, PIN, QIN) RESULT (RES)
!***BEGIN PROLOGUE  DBETAI
!***PURPOSE  Calculate the incomplete Beta function.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7F
!***TYPE      DOUBLE PRECISION (BETAI-S, DBETAI-D)
!***KEYWORDS  FNLIB, INCOMPLETE BETA FUNCTION, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
!   DBETAI calculates the DOUBLE PRECISION incomplete beta function.
!
!   The incomplete beta function ratio is the probability that a
!   random variable from a beta distribution having parameters PIN and
!   QIN will be less than or equal to X.
!
!     -- Input Arguments -- All arguments are DOUBLE PRECISION.
!   X      upper limit of integration.  X must be in (0,1) inclusive.
!   PIN    first beta distribution parameter.  PIN must be .GT. 0.0.
!   QIN    second beta distribution parameter.  QIN must be .GT. 0.0.
!
!***REFERENCES  Nancy E. Bosten and E. L. Battiste, Remark on Algorithm
!                 179, Communications of the ACM 17, 3 (March 1974),
!                 pp. 156.
!***ROUTINES CALLED  D1MACH, DLBETA, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   890911  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920528  DESCRIPTION and REFERENCES sections revised.  (WRB)
!***END PROLOGUE  DBETAI
!      DOUBLE PRECISION X, PIN, QIN, ALNEPS, ALNSML, C, EPS, FINSUM, P, &
!                       PS, Q, SML, TERM, XB, XI, Y, DLBETA, P1
!      LOGICAL FIRST
!      INTEGER I,IB,N
!      DOUBLE PRECISION RES
!      SAVE EPS, ALNEPS, SML, ALNSML, FIRST
!      DATA FIRST /.TRUE./
!!***FIRST EXECUTABLE STATEMENT  DBETAI
!      IF (FIRST) THEN
!         EPS = D1MACH(3)
!         ALNEPS = LOG (EPS)
!         SML = D1MACH(1)
!         ALNSML = LOG (SML)
!      ENDIF
!      FIRST = .FALSE.

!      IF (X .LT. 0.D0 .OR. X .GT. 1.D0) THEN
!         CALL STOP_APPLICATION('DEBATI : X IS NOT IN THE RANGE (0,1)')
!      END IF
!      IF (PIN .LE. 0.D0 .OR. QIN .LE. 0.D0) THEN
!        CALL STOP_APPLICATION ('DBETAI : P AND/OR Q IS LE ZERO')
!      END IF

!      Y = X
!      P = PIN
!      Q = QIN
!      IF (Q.LE.P .AND. X.LT.0.8D0) GO TO 20
!      IF (X.LT.0.2D0) GO TO 20
!      Y = 1.0D0 - Y
!      P = QIN
!      Q = PIN

! 20   IF ((P+Q)*Y/(P+1.D0).LT.EPS) GO TO 80

!! EVALUATE THE INFINITE SUM FIRST.  TERM WILL EQUAL
!! Y**P/BETA(PS,P) * (1.-PS)-SUB-I * Y**I / FAC(I) .

!      PS = Q - AINT(Q)
!      IF (PS.EQ.0.D0) PS = 1.0D0
!      XB = P*LOG(Y) - DLBETA(PS,P) - LOG(P)
!      RES = 0.0D0
!      IF (XB.LT.ALNSML) GO TO 40

!      RES = EXP (XB)
!      TERM = RES*P
!      IF (PS.EQ.1.0D0) GO TO 40
!      N = MAX (ALNEPS/LOG(Y), 4.0D0)
!      DO 30 I=1,N
!        XI = I
!        TERM = TERM * (XI-PS)*Y/XI
!        RES = RES + TERM/(P+XI)
! 30   CONTINUE

!! NOW EVALUATE THE FINITE SUM, MAYBE.

! 40   IF (Q.LE.1.0D0) GO TO 70

!      XB = P*LOG(Y) + Q*LOG(1.0D0-Y) - DLBETA(P,Q) - LOG(Q)
!      IB = MAX (XB/ALNSML, 0.0D0)
!      TERM = EXP(XB - IB*ALNSML)
!      C = 1.0D0/(1.D0-Y)
!      P1 = Q*C/(P+Q-1.D0)

!      FINSUM = 0.0D0
!      N = Q
!      IF (Q.EQ.DBLE(N)) N = N - 1
!      DO 50 I=1,N
!        IF (P1.LE.1.0D0 .AND. TERM/EPS.LE.FINSUM) GO TO 60
!        XI = I
!        TERM = (Q-XI+1.0D0)*C*TERM/(P+Q-XI)

!        IF (TERM.GT.1.0D0) IB = IB - 1
!        IF (TERM.GT.1.0D0) TERM = TERM*SML

!        IF (IB.EQ.0) FINSUM = FINSUM + TERM
! 50   CONTINUE

! 60   RES = RES + FINSUM
! 70   IF (Y.NE.X .OR. P.NE.PIN) RES = 1.0D0 - RES
!      RES = MAX (MIN (RES, 1.0D0), 0.0D0)
!      RETURN

! 80   RES = 0.0D0
!      XB = P*LOG(MAX(Y,SML)) - LOG(P) - DLBETA(P,Q)
!      IF (XB.GT.ALNSML .AND. Y.NE.0.0D0) RES = EXP(XB)
!      IF (Y.NE.X .OR. P.NE.PIN) RES = 1.0D0 - RES

!      RETURN
!      END FUNCTION DBETAI

      subroutine ginv1(a,n,m,tol,irank)
! returns generalized inverse of matrix x of size n x n declared
! as m x m. tol is working zero and irank returns the rank of
! the matrix. w is a work vector of size m,
! by rohan fernando, slightly structured by i. misztal 05/05/87

      integer :: n,m
      double precision a(m,m),w(m),re,sum,tol
      integer :: i,j,iim1,ii,iii,irank
      irank=n
      do 10 i=1,n
         do 20 j=1,i-1
              re=a(i,j)
              do 20 ii=i,n
20                 a(ii,i)=a(ii,i)-re*a(ii,j)
         if (a(i,i).lt.tol) then
              a(i,i)=0.0
              do 45 ii=i+1,n
45                 a(ii,i)=0.
           irank=irank-1
           else
              a(i,i)=sqrt(a(i,i))
              do 40 ii=i+1,n
40                a(ii,i)=a(ii,i)/a(i,i)
         endif
10    continue

      do 100 i=1,n
         if (a(i,i).eq.0.) then
              do 150 ii=i+1,n
150                a(ii,i)=0.
           else
              a(i,i)=1.0/ a(i,i)
              do 200 ii=i+1,n
200               w(ii)=0.0
              do 300 ii=i+1,n
                  iim1=ii-1
                  re=a(iim1,i)
                  do 400 iii=ii,n
400                   w(iii)=w(iii)-a(iii,iim1)*re
                  if (a(ii,ii).eq.0.) then
                      a(ii,i)=0.
                    else
                      a(ii,i)=w(ii)/a(ii,ii)
                  endif
300           continue
          endif
100     continue

      do 110 j=1,n
         do 110 i=j,n
              sum=0
              do 130 ii=i,n
130                sum=sum+a(ii,j)*a(ii,i)
110           a(i,j)=sum
      do 600 i=1,n
          do 600 j=i,n
600           a(i,j)=a(j,i)
      return
      end subroutine ginv1


      SUBROUTINE tdiag (n, tol, a, d, e, z, maxdim, ifault)

! Code converted using TO_F90 by Alan Miller
! Date: 2003-12-03  Time: 12:16:20

INTEGER, INTENT(IN)      :: n
REAL (dp), INTENT(OUT)   :: tol
INTEGER, INTENT(IN)  :: maxdim
REAL (dp), INTENT(IN)    :: a(maxdim,maxdim)
REAL (dp), INTENT(OUT)   :: d(maxdim)
REAL (dp), INTENT(OUT)   :: e(maxdim)
REAL (dp), INTENT(OUT)   :: z(maxdim,maxdim)
INTEGER, INTENT(OUT)     :: ifault

REAL (dp), PARAMETER  :: eta = 1.0D-37, eps = 1.0D-14

!    Algorithm as 60.1 appl.statist. (1973) vol.22 no.2

!    Reduces real symmetric matrix to tridiagonal form

!    tol is a machine dependent constant , tol = eta/eps , where
!    eta is the smallest positive number representable in the
!    computer and eps is the smallest positive number for which
!    1+eps.ne.1.

!      eta=eps*tol
!      eps=0.7105427358e-14
!      tol=0.3131513063e-293
!      precis=1.0e-14

!      NB
!        Real constants must be <= 15 decimal digits
!        the range of a real constant is from 1.0e-293 to 1.0e+322

INTEGER    :: i, i1, j, j1, k, l
REAL (dp)  :: f, g, h, hh,zero,one

zero=0.d0
one=1.d0
tol=eta/eps
ifault=1
IF (n <= 1) RETURN
ifault=0
DO  i=1,n
  z(i,1:i)=a(i,1:i)
END DO
i=n
DO  i1=2,n
  l=i-2
  f=z(i,i-1)
  g=zero
  IF (l < 1) GO TO 30
  g = SUM( z(i,1:l)**2 )
  30 h=g + f*f

!       If g is too small for orthogonality to be guaranteed, the
!       transformation is skipped

  IF (g > tol) GO TO 40
  e(i)=f
  d(i)=zero
  GO TO 100

  40 l=l+1
  g=SQRT(h)
  IF (f >= zero) g=-g
  e(i)=g
  h=h - f*g
  z(i,i-1)=f-g
  f=zero
  DO  j=1,l
    z(j,i)=z(i,j)/h
    g=zero

!       Form element of a * u

    DO  k=1,j
      g=g + z(j,k)*z(i,k)
    END DO
    IF (j >= l) GO TO 70
    j1=j+1
    DO  k=j1,l
      g=g + z(k,j)*z(i,k)
    END DO

!       Form element of p

    70 e(j)=g/h
    f=f + g*z(j,i)
  END DO

!       Form k

  hh=f/(h+h)

!       Form reduced a

  DO  j=1,l
    f=z(i,j)
    g=e(j)-hh*f
    e(j)=g
    DO  k=1,j
      z(j,k)=z(j,k) - f*e(k) - g*z(i,k)
    END DO
  END DO
  d(i)=h
  100 i=i-1
END DO
d(1)=zero
e(1)=zero

!       Accumulation of transformation matrices

DO  i=1,n
  l=i-1
  IF (d(i) == zero .OR. l == 0) GO TO 140
  DO  j=1,l
    g=zero
    DO  k=1,l
      g=g + z(i,k)*z(k,j)
    END DO
    DO  k=1,l
      z(k,j)=z(k,j) - g*z(k,i)
    END DO
  END DO
  140 d(i)=z(i,i)
  z(i,i)=one
  IF (l == 0) CYCLE
  DO  j=1,l
    z(i,j)=zero
    z(j,i)=zero
  END DO
END DO
RETURN
END SUBROUTINE tdiag


SUBROUTINE lrvt (n, precis, d, e, z, ifault, maxdim)

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(OUT)     :: precis
INTEGER, INTENT(IN)        :: maxdim
REAL (dp), INTENT(IN OUT)  :: d(maxdim)
REAL (dp), INTENT(OUT)     :: e(maxdim)
REAL (dp), INTENT(IN OUT)  :: z(maxdim,maxdim)
INTEGER, INTENT(OUT)       :: ifault

!    Algorithm AS 60.2 appl.statist. (1973) vol.22, no.2

!    Finds latent roots and vectors of tridiagonal matrix

INTEGER, PARAMETER    :: mits = 30
REAL (dp), PARAMETER  :: two = 2.0D0

INTEGER    :: i, i1, j, jj, k, l, m, m1, n1
REAL (dp)  :: b, c, f, g, h, p, pr, r, s,zero,one
zero=0.d0
one=1.d0
precis=1.0D-14
ifault=2
IF (n <= 1) RETURN
ifault=1
n1=n-1
DO  i=2,n
  e(i-1)=e(i)
END DO
e(n)=zero
b=zero
f=zero
DO  l=1,n
  jj=0
  h=precis*(ABS(d(l))+ABS(e(l)))
  IF (b < h) b=h

!       Look for small sub-diagonal element

  DO  m1=l,n
    m=m1
    IF (ABS(e(m)) <= b) EXIT
  END DO
  IF (m == l) GO TO 90
  40 IF (jj == mits) RETURN
  jj=jj+1

!       Form shift

  p=(d(l+1)-d(l))/(two*e(l))
  r=SQRT(p*p + one)
  pr=p+r
  IF (p < zero) pr=p-r
  h=d(l)-e(l)/pr
  DO  i=l,n
    d(i)=d(i)-h
  END DO
  f=f+h

!       QL transformation

  p=d(m)
  c=one
  s=zero
  m1=m-1
  i=m
  DO  i1=l,m1
    j=i
    i=i-1
    g=c*e(i)
    h=c*p
    IF (ABS(p) >= ABS(e(i))) GO TO 60
    c=p/e(i)
    r=SQRT(c*c + one)
    e(j)=s*e(i)*r
    s=one/r
    c=c/r
    GO TO 70
    60 c=e(i)/p
    r=SQRT(c*c + one)
    e(j)=s*p*r
    s=c/r
    c=one/r
    70 p=c*d(i) - s*g
    d(j)=h + s*(c*g + s*d(i))

!       Form vector

    DO  k=1,n
      h=z(k,j)
      z(k,j)=s*z(k,i) + c*h
      z(k,i)=c*z(k,i) - s*h
    END DO
  END DO
  e(l)=s*p
  d(l)=c*p
  IF (ABS(e(l)) > b) GO TO 40
  90 d(l)=d(l) + f
END DO

!       Order latent roots and vectors

DO  i=1,n1
  k=i
  p=d(i)
  i1=i+1
  DO  j=i1,n
    IF (d(j) <= p) CYCLE
    k=j
    p=d(j)
  END DO
  IF (k == i) CYCLE
  d(k)=d(i)
  d(i)=p
  DO  j=1,n
    p=z(j,i)
    z(j,i)=z(j,k)
    z(j,k)=p
  END DO
END DO
ifault=0
RETURN
END SUBROUTINE lrvt


end module m_qtlmap_math




