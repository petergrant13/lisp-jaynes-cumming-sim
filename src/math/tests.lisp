;; This file has a bunch of unit tests for the matrix operations code

(defun assert-near (a b &optional (eps 1e-8))
  (unless (< (abs (- a b)) eps)
    (error "Assertion failed: ~A != ~A" a b)))

(defmacro assert-signals (condition-type &body body)
  `(handler-case
       (progn
         ,@body
         (error "Expected condition ~A was not signaled"
                ',condition-type))

     (,condition-type (c)
       t)))

(defun matrix-equal-near (A B &optional (eps 1e-12))
  (let ((rows (array-dimension A 0))
        (cols (array-dimension A 1)))
    (unless (and (= rows (array-dimension B 0))
                 (= cols (array-dimension B 1)))
      (return-from matrix-equal-near nil))
    
    (dotimes (i rows t)
      (dotimes (j cols)
        (unless (< (abs (- (aref A i j)
                           (aref B i j)))
                   eps)
          (return-from matrix-equal-near nil))))))

(defun test-gauss-jordan-identity ()
  (let ((A (identity-matrix 3)))
    (gauss-jordan A)
    (unless (matrix-equal-near A (identity-matrix 3))
      (error "Identity matrix test failed")))
  
  (format t "Identity test passed.~%"))

(defun test-solve-system ()
  (let ((A (make-matrix 2 2 0.0))
        (b (make-vector 2 0.0)))

    ;; A matrix
    (setf (aref A 0 0) 1.0)
    (setf (aref A 0 1) 1.0)
    (setf (aref A 1 0) 1.0)
    (setf (aref A 1 1) -1.0)

    ;; RHS vector
    (setf (aref b 0) 2.0)
    (setf (aref b 1) 0.0)

    (let ((x (solve-system A b)))
      (assert-near (aref x 0 0) 1.0)
      (assert-near (aref x 1 0) 1.0)))

  (format t "Linear solve test passed.~%"))

(defun test-singular-system ()
  (let ((A (make-matrix 2 2 0.0))
        (b (make-vector 2 0.0)))

    (setf (aref A 0 0) 1.0)
    (setf (aref A 0 1) 2.0)
    (setf (aref A 1 0) 2.0)
    (setf (aref A 1 1) 4.0)

    (assert-signals error
      (solve-system A b))

    (format t "Singular matrix test passed.~%")))

(defun run-tests ()
  (test-gauss-jordan-identity)
  (test-solve-system)
  (test-singular-system)

  (format t "~%All tests passed.~%"))
