;These are helper functions for complex math :

(defun absolute-square (z)
  (expt (abs z) 2))
  
(defun phase-degrees (z)
  (* 180 (/ (phase z) pi)))
