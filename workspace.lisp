;;; load program
(load "/home/wiz/ngnet-online-em/ngnet-online-em.lisp")
(asdf:oos 'asdf:load-op :cl-store)

;;; データ点
(defparameter x_in 1d0)
(defparameter x_out 3d0)

;;; データセットを生成
(defparameter x_in.dat (mapcar (lambda (x) (random-uniform (- (* 1.5 pi)) (* 1.5 pi))) #[0 20000]))
(defparameter x_out.dat (mapcar (lambda (x) (+ (sin x) (random-normal :sd 0.3))) x_in.dat))

;; (cl-store:store x_in.dat "/home/wiz/research/NGnet/oibem/data/x_in.1dim.data1")
;; (cl-store:store x_in.dat "/home/wiz/research/NGnet/oibem/data/x_in.1dim.data2")
;; (cl-store:store x_out.dat "/home/wiz/research/NGnet/oibem/data/x_out.1dim.data1")
;; (cl-store:store x_out.dat "/home/wiz/research/NGnet/oibem/data/x_out.1dim.data2")

(defvar x_in.dat1)
(defvar x_in.dat2)
(defvar x_out.dat1)
(defvar x_out.dat2)
(setf x_in.dat1 (cl-store:restore "/home/wiz/research/NGnet/oibem/data/x_in.1dim.data1"))
(setf x_in.dat2 (cl-store:restore "/home/wiz/research/NGnet/oibem/data/x_in.1dim.data2"))
(setf x_out.dat1 (cl-store:restore "/home/wiz/research/NGnet/oibem/data/x_out.1dim.data1"))
(setf x_out.dat2 (cl-store:restore "/home/wiz/research/NGnet/oibem/data/x_out.1dim.data2"))

(defvar data1)
(defvar data2)
(setf data1 (make-data-stream x_in.dat1 x_out.dat1))
(setf data2 (make-data-stream x_in.dat2 x_out.dat2))

;;; データセットをプロット
(plot-lst x_out.dat :x-lst x_in.dat :style nil :title "data-set" :xlabel "x_{in}" :ylabel "x_{out}")
(plot-lst x_out.dat1 :x-lst x_in.dat1 :style nil :title "data-set" :xlabel "x_{in}" :ylabel "x_{out}")
(plot-lst x_out.dat2 :x-lst  x_in.dat2 :style nil :title "data-set" :xlabel "x_in" :ylabel "x_out")

(defparameter data-stream
  (make-instance 'data-stream-general-source
		 :input-generator
		 (lambda () (random-uniform (- (* 1.5 pi)) (* 1.5 pi)))
		 :input-output-function
		 (lambda (x_in) (+ (sin x_in) (random-normal :sd 0.3)))))

;;; 同時分布をプロット
(defun plot-joint-1dim (session xrange yrange)
  (plot-lst-2dim
   (lambda (x_in x_out)
     (ssum (apply #'joint-probability x_in x_out (session-params session))))
   xrange yrange :surface nil :xlabel "x_{in}" :ylabel "x_{out}" :title "P(x_{in},x_{out})"
   ))

;;; 事後分布をプロット
(defun plot-post-1dim (session xrange yrange)
  (plot-lst-2dim
   (lambda (x_in x_out)
     (let ((joint (ssum (apply #'joint-probability x_in x_out (session-params session)))))
       (/ joint (sum-denomi x_in yrange session))))
   xrange yrange :surface nil :xlabel "x_{in}" :ylabel "x_{out}" :title "P(x_{out}|x_{in})"
   ))

(defun sum-denomi (x y.dat session)
  (ssum (mapcar (lambda (x_out)
		  (ssum (apply #'joint-probability x x_out (session-params session))))
		y.dat)))

;;; 予測線(初期パラメータ)をプロット
(defun plot-prediction-1dim (session xrange)
  (plot-lst (mapcar (lambda (x) (prediction x session))
		    xrange)))

;; sessionの定義とNGnetの学習
(defparameter session (make-initialized-session 100 1 1 :lambda-factor 0.999d0 :alpha 0.1d0 :LL 0d0))
(defparameter learned-session (oem-data-list x_in.dat x_out.dat session))

(defparameter session2 (oem-1step (car x_in.dat) (car x_out.dat) session))
(defparameter session3 (oem-1step (cadr x_in.dat) (cadr x_out.dat) session2))

(joint-probability (cadr x_in.dat) (cadr x_out.dat) (session-mu session2) (session-sigma session2) '(1 1 1) (session-W session2) (session-s session2))

(oem-data-stream data-stream session)
(plot-prediction-1dim learned-session (seq (- (* 1.5 pi)) (* 1.5 pi) :by 0.05)) ;*last-session*

;; 入出力の同時分布
(plot-joint-1dim *last-session* ;learned-session
		 (seq (- (* 1.5 pi)) (* 1.5 pi) :by 0.1)
		 (seq (- 2.0) 2.0 :by 0.1))


;; predict-samplingの結果をプロット(正弦曲線の周りにうまく散らばっているか)
(let ((x-lst (seq (- (* 2.0 pi)) (* 2.0 pi) :by 0.01)))
  (plot-lst (mapcar (lambda (x_in)
		      (predict-sampling x_in learned-session))
		    x-lst)
	    :x-lst x-lst :style "points" :yrange '(-2.0 2.0)))

(let ((x-lst (seq (- (* 2.0 pi)) (* 2.0 pi) :by 0.01)))
  (plot-lst (mapcar (lambda (x_in)
		      (+ (prediction x_in learned-session)
			 (random-normal :sd 0.5)))
		    x-lst)
	    :x-lst x-lst :style "points"))