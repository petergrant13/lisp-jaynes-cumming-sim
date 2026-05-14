(defpackage :matrix-operations
  (:use :cl)
  (:export
   *eps*
   ;; constructors
   make-matrix make-vector
   print-matrix identity-matrix copy-matrix

   ;; vector operations
   vector-add vector-sub vector-scale
   dot-product vector-norm

   ;; matrix arithmetic
   matrix-add matrix-sub matrix-scale
   matrix-multiply matrix-vector-multiply
   transpose trace

   ;; row operations
   swap-rows scale-row add-rows

   ;; linear solver
   gauss-jordan solve-system))

(in-package :matrix-operations)

(defparameter *eps* 1e-12)

;;; --------------------------
;;; Constructors
;;; --------------------------

(defun make-matrix (rows cols &optional (init 0))
  (make-array (list rows cols) :initial-element init))
;; update element with (setf (aref M i j) value)
;; M = matrix
;; i = row index
;; j = column index
;; Indexing starts at 0.


(defun make-vector (n &optional (init 0))
  (make-array n :initial-element init))

(defun copy-matrix (A)
  (let* ((rows (array-dimension A 0))
         (cols (array-dimension A 1))
         (B (make-matrix rows cols)))
    (dotimes (i rows)
      (dotimes (j cols)
        (setf (aref B i j) (aref A i j))))
    B))

(defun print-matrix (M)
  (dotimes (i (array-dimension M 0))
    (dotimes (j (array-dimension M 1))
      (format t "~8,4f " (float (aref M i j))))
    (terpri)))

(defun identity-matrix (n)
  (let ((m (make-matrix n n)))
    (dotimes (i n)
      (setf (aref m i i) 1))
    m))

;;; --------------------------
;;; Vector operations
;;; --------------------------

(defun vector-add (a b)
  (unless (= (length a) (length b))
    (error "Vector dimensions do not match"))
  (let ((v (make-vector (length a))))
    (dotimes (i (length a))
      (setf (aref v i) (+ (aref a i) (aref b i))))
    v))

(defun vector-sub (a b)
  (unless (= (length a) (length b))
    (error "Vector dimensions do not match"))
  (let ((v (make-vector (length a))))
    (dotimes (i (length a))
      (setf (aref v i) (- (aref a i) (aref b i))))
    v))

(defun vector-scale (s v)
  (let ((out (make-vector (length v))))
    (dotimes (i (length v))
      (setf (aref out i) (* s (aref v i))))
    out))

(defun dot-product (a b)
  (unless (= (length a) (length b))
    (error "Vector dimensions do not match"))
  (let ((sum 0))
    (dotimes (i (length a))
      (incf sum (* (aref a i) (aref b i))))
    sum))

(defun vector-norm (v)
  (sqrt (dot-product v v)))

;;; --------------------------
;;; Matrix arithmetic
;;; --------------------------

(defun matrix-add (A B)
  (unless (and (= (array-dimension A 0) (array-dimension B 0))
               (= (array-dimension A 1) (array-dimension B 1)))
    (error "Matrix dimensions do not match"))
  (let* ((rows (array-dimension A 0))
         (cols (array-dimension A 1))
         (C (make-matrix rows cols)))
    (dotimes (i rows)
      (dotimes (j cols)
        (setf (aref C i j) (+ (aref A i j) (aref B i j)))))
    C))

(defun matrix-sub (A B)
  (unless (and (= (array-dimension A 0) (array-dimension B 0))
               (= (array-dimension A 1) (array-dimension B 1)))
    (error "Matrix dimensions do not match"))
  (let* ((rows (array-dimension A 0))
         (cols (array-dimension A 1))
         (C (make-matrix rows cols)))
    (dotimes (i rows)
      (dotimes (j cols)
        (setf (aref C i j) (- (aref A i j) (aref B i j)))))
    C))

(defun matrix-scale (s A)
  (let* ((rows (array-dimension A 0))
         (cols (array-dimension A 1))
         (C (make-matrix rows cols)))
    (dotimes (i rows)
      (dotimes (j cols)
        (setf (aref C i j) (* s (aref A i j)))))
    C))

(defun matrix-multiply (A B)
  (unless (= (array-dimension A 1) (array-dimension B 0))
    (error "Matrix dimensions incompatible"))
  (let* ((rows (array-dimension A 0))
         (cols (array-dimension B 1))
         (inner (array-dimension A 1))
         (C (make-matrix rows cols 0)))
    (dotimes (i rows)
      (dotimes (j cols)
        (dotimes (k inner)
          (incf (aref C i j)
                (* (aref A i k) (aref B k j))))))
    C))

(defun matrix-vector-multiply (A v)
  (unless (= (array-dimension A 1) (length v))
    (error "Matrix/vector dimensions incompatible"))
  (let* ((rows (array-dimension A 0))
         (result (make-vector rows)))
    (dotimes (i rows)
      (let ((sum 0))
        (dotimes (j (length v))
          (incf sum (* (aref A i j) (aref v j))))
        (setf (aref result i) sum)))
    result))

(defun transpose (M)
  (let* ((rows (array-dimension M 0))
         (cols (array-dimension M 1))
         (MT (make-matrix cols rows)))
    (dotimes (i rows)
      (dotimes (j cols)
        (setf (aref MT j i) (aref M i j))))
    MT))

(defun trace (M)
  (unless (= (array-dimension M 0) (array-dimension M 1))
    (error "Trace requires square matrix"))
  (let ((sum 0))
    (dotimes (i (array-dimension M 0))
      (incf sum (aref M i i)))
    sum))

;;; --------------------------
;;; Row operations
;;; --------------------------

(defun swap-rows (M i j)
  (let ((cols (array-dimension M 1)))
    (dotimes (k cols)
      (rotatef (aref M i k) (aref M j k)))))

(defun scale-row (M i factor)
  (let ((cols (array-dimension M 1)))
    (dotimes (k cols)
      (setf (aref M i k) (* factor (aref M i k))))))

(defun add-rows (M src dest factor)
  (let ((cols (array-dimension M 1)))
    (dotimes (k cols)
      (incf (aref M dest k) (* factor (aref M src k))))))

;;; --------------------------
;;; Gauss–Jordan elimination
;;; --------------------------

(defun gauss-jordan (M)
  (let* ((rows (array-dimension M 0))
         (cols (array-dimension M 1))
         (lead 0)
         (rank 0))

    (dotimes (r rows)

      ;; stop if we've exhausted columns
      (when (>= lead cols)
        (return))

      ;; find pivot row
      (let ((pivot-row nil)
            (max-val 0.0))

        (dotimes (k (- rows r))
          (let* ((ii (+ r k))
                 (val (abs (aref M ii lead))))
            (when (> val max-val)
              (setf max-val val
                    pivot-row ii))))

        ;; if no pivot exists in this column,
        ;; move to next column and retry same row
        (if (or (null pivot-row)
                (<= max-val *eps*))

            (progn
              (incf lead)
              (decf r))

            (progn

              ;; swap pivot into place
              (when (/= pivot-row r)
                (swap-rows M r pivot-row))

              ;; normalize pivot row
              (let ((pivot (aref M r lead)))
                (when (<= (abs pivot) *eps*)
                  (error "Singular matrix"))

                (scale-row M r (/ 1.0 pivot)))

              ;; eliminate all other rows
              (dotimes (j rows)
                (unless (= j r)
                  (let ((factor (aref M j lead)))
                    (when (> (abs factor) *eps*)
                      (add-rows M r j (- factor))))))

              (incf rank)
              (incf lead)))))

    ;; detect inconsistent rows
    ;;
    ;; [0 0 0 | nonzero]
    ;;
    ;; means no solution exists
    (dotimes (i rows)

      (let ((all-zero t))

        ;; check coefficient part only
        (dotimes (j (1- cols))
          (when (> (abs (aref M i j)) *eps*)
            (setf all-zero nil)))

        ;; inconsistent row
        (when (and all-zero
                   (> (abs (aref M i (1- cols))) *eps*))
          (error "Inconsistent system: no solution"))))

    ;; detect non-unique solutions
    ;;
    ;; rank < number of variables
    ;;
    (when (< rank (1- cols))
      (error "Singular or underdetermined system"))

    M))

(defun solve-system (A b)
  (let* ((rows (array-dimension A 0))
         (cols (array-dimension A 1)))
    (let* ((M (make-matrix rows (1+ cols))))
      (dotimes (i rows)
        (dotimes (j cols)
          (setf (aref M i j) (aref A i j)))
        (setf (aref M i cols)
              (if (= (array-rank b) 1)
                  (aref b i)
                  (aref b i 0))))
      (gauss-jordan M)
      (let ((x (make-matrix cols 1)))
        (dotimes (i cols)
          (setf (aref x i 0) (aref M i cols)))
        x))))
