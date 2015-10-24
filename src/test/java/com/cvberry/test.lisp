(defun element-wise-matrix (fn A B)
  (let* ((len (array-total-size A))
         (m   (car (array-dimensions A)))
         (n   (cadr (array-dimensions A)))
         (C   (make-array `(,m ,n) :initial-element 0.0d0)))
 
    (loop for i from 0 to (1- len) do
         (setf (row-major-aref C i) 
               (funcall fn
                        (row-major-aref A i)
                        (row-major-aref B i))))
    C))
 
;; A.+B, A.-B, A.*B, A./B, A.^B.
(defun m+ (A B) (element-wise-matrix #'+    A B))
(defun m- (A B) (element-wise-matrix #'-    A B))
(defun m* (A B) (element-wise-matrix #'*    A B))
(defun m/ (A B) (element-wise-matrix #'/    A B))
(defun m^ (A B) (element-wise-matrix #'expt A B))

(defun element-wise-scalar (fn A c)
  (let* ((len (array-total-size A))
         (m   (car (array-dimensions A)))
         (n   (cadr (array-dimensions A)))
         (B   (make-array `(,m ,n) :initial-element 0.0d0)))
 
    (loop for i from 0 to (1- len) do
         (setf (row-major-aref B i) 
               (funcall fn
                        (row-major-aref A i)
                        c)))
    B))
 
;; c.+A, A.-c, c.*A, A./c, A.^c.
(defun .+ (c A) (element-wise-scalar #'+    A c))
(defun .- (A c) (element-wise-scalar #'-    A c))
(defun .* (c A) (element-wise-scalar #'*    A c))
(defun ./ (A c) (element-wise-scalar #'/    A c))
(defun .^ (A c) (element-wise-scalar #'expt A c))

(defun mmul (A B)
  (let* ((m (car (array-dimensions A)))
         (n (cadr (array-dimensions A)))
         (l (cadr (array-dimensions B)))
         (C (make-array `(,m ,l) :initial-element 0)))
    (loop for i from 0 to (- m 1) do
              (loop for k from 0 to (- l 1) do
                    (setf (aref C i k)
                          (loop for j from 0 to (- n 1)
                                sum (* (aref A i j)
                                       (aref B j k))))))
    C))

;; Transpose a mxn matrix A to a nxm matrix B=A'.
(defun mtp (A)
  (let* ((m (array-dimension A 0))
         (n (array-dimension A 1))
         (B (make-array `(,n ,m) :initial-element 0)))
    (loop for i from 0 below m do
          (loop for j from 0 below n do
                (setf (aref B j i)
                      (aref A i j))))
    B))


(defun sign (x)
  (if (zerop x)
      x
      (/ x (abs x))))

(defun norm (x)
  (let ((len (car (array-dimensions x))))
    (sqrt (loop for i from 0 to (1- len) sum (expt (aref x i 0) 2)))))

(norm #2a((1 2 3)
                   (4 5 400)))


(elt #1a(1 2 3 4) 2)
 
(defun make-unit-vector (dim)
  (let ((vec (make-array `(,dim ,1) :initial-element 0.0d0)))
    (setf (aref vec 0 0) 1.0d0)
    vec))
 
;; Return a nxn identity matrix.
(defun eye (n)
  (let ((I (make-array `(,n ,n) :initial-element 0)))
    (loop for j from 0 to (- n 1) do
          (setf (aref I j j) 1))
    I))
 
(defun array-range (A ma mb na nb)
  (let* ((mm (1+ (- mb ma)))
         (nn (1+ (- nb na)))
         (B (make-array `(,mm ,nn) :initial-element 0.0d0)))
 
    (loop for i from 0 to (1- mm) do
         (loop for j from 0 to (1- nn) do
              (setf (aref B i j)
                    (aref A (+ ma i) (+ na j)))))
    B))
 
(defun rows (A) (car  (array-dimensions A)))
(defun cols (A) (cadr (array-dimensions A)))
(defun mcol (A n) (array-range A 0 (1- (rows A)) n n))
(defun mrow (A n) (array-range A n n 0 (1- (cols A))))
 
(defun array-embed (A B row col)
  (let* ((ma (rows A))
         (na (cols A))
         (mb (rows B))
         (nb (cols B))
         (C  (make-array `(,ma ,na) :initial-element 0.0d0)))
 
    (loop for i from 0 to (1- ma) do
         (loop for j from 0 to (1- na) do
              (setf (aref C i j) (aref A i j))))
 
    (loop for i from 0 to (1- mb) do
         (loop for j from 0 to (1- nb) do
              (setf (aref C (+ row i) (+ col j))
                    (aref B i j))))
 
    C))
 
(defun make-householder (a)
  (let* ((m    (car (array-dimensions a)))
         (s    (sign (aref a 0 0)))
         (e    (make-unit-vector m))
         (u    (m+ a (.* (* (norm a) s) e)))
         (v    (./ u (aref u 0 0)))
         (beta (/ 2 (aref (mmul (mtp v) v) 0 0))))
 
    (m- (eye m)
        (.* beta (mmul v (mtp v))))))
 
(defun qr (A)
  (let* ((m (car  (array-dimensions A)))
         (n (cadr (array-dimensions A)))
         (Q (eye m)))
 
    ;; Work on n columns of A.
    (loop for i from 0 to (if (= m n) (- n 2) (- n 1)) do
 
         ;; Select the i-th submatrix. For i=0 this means the original matrix A.
         (let* ((B (array-range A i (1- m) i (1- n)))
                ;; Take the first column of the current submatrix B.
                (x (mcol B 0))
                ;; Create the Householder matrix for the column and embed it into an mxm identity.
                (H (array-embed (eye m) (make-householder x) i i)))
 
           ;; The product of all H matrices from the right hand side is the orthogonal matrix Q.
           (setf Q (mmul Q H))
 
           ;; The product of all H matrices with A from the LHS is the upper triangular matrix R.
           (setf A (mmul H A))))
 
    ;; Return Q and R.
    (values Q A)))

(qr #2A((12 -51 4) (6 167 -68) (-4 24 -41)))
