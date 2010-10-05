;;; 入出力の同時分布
(plot-lst-2dim
 (lambda (x y) (ssum (apply #'joint-probability x y (session-params session))))
 (seq -4.0d0 4.0d0 :by 0.1) (seq -4.0d0 4.0d0 :by 0.1) :surface nil)

;;; 入力が与えられたときの出力の予測分布
(plot-lst-2dim
 (lambda (x y)
   (/ (ssum
       (mapcar (lambda (unit)
		 (gaussian y
			   (aref (m* (nth unit (session-W session)) (ext-vector x)) 0 0)
			   (nth unit (session-s session))))
	       #[0 2]))
      3))
 (seq -4.0d0 4.0d0 :by 0.1) (seq -4.0d0 4.0d0 :by 0.1) :surface nil)

;;; 予測線
(plot-lst (mapcar (lambda (x) (prediction x session)) (seq -4.0d0 4.0d0 :by 0.1d0)))

;;; 予測分布の断面
(let ((x -0.2d0))
  (plot-lst
   (mapcar
    (lambda (y)
      (/ (ssum
	  (mapcar (lambda (unit)
		    (gaussian y
			      (aref (m* (nth unit (session-W session)) (ext-vector x)) 0 0)
			      (nth unit (session-s session))))
		  #[0 2]))
	 3))
    (seq -4.0d0 4.0d0 :by 0.1))))

;;; 予測分布からのサンプリング
(setf predict-points (scollect 10000 (lambda () (predict-sampling -2d0 session))))
(histogram predict-points 20)

