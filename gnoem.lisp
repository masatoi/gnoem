;;; -*- Coding: utf-8; Mode: Lisp -*-

;;; (proclaim '(optimize speed))
(declaim (optimize speed))

(in-package :gnoem)
;; ;;; load Gaussian functions
;; (load "/home/wiz/program/cl/multivariate-gaussian.lisp")

;; ;;; load NLISP and define plot function
;; (load "/home/wiz/program/cl/wiz-util/frontend-nlisp.lisp")


;;; NGnet =============================================================

;; make extend vector x~ = [x, 1]^T
(defun ext-vector (x)
  (if (arrayp x)
      (vector-cat x (make-vector 1 :initial-element 1.0d0))
      (make-vector 2 :initial-contents (list x 1.0d0))))

;;; joint distribution of NGnet ; P(x, y, i | parameter)
(defun joint-probability (x y mu Sigma Sigma-inv W s <<1>>-acc Lambda~-acc)
  (let ((x~ (ext-vector x))
	(|t| (length mu))) ; number of units
    (mapcar
     (lambda (mu Sigma Sigma-inv W s <<1>> Lambda~)
       (setf *last-x* x
	     *last-y* y
	     *last-mu* mu
	     *last-Sigma* Sigma
	     *last-Sigma-inv* Sigma-inv
	     *last-W* W
	     *last-s* s
	     *last-<<1>>* <<1>>
	     *last-Lambda~* Lambda~)
       (let ((prob-val
	      (/ (* (universal-gaussian x mu Sigma) ; :inv-Sigma Sigma-inv
		    (if (arrayp y) ; 出力が多次元かどうかで 場合分け
			(let ((out-sd (m* s (umat (array-dimension y 0)))))
			  (universal-gaussian y (m* W x~) out-sd :inv-Sigma out-sd))
			(universal-gaussian y (aref (m* W x~) 0 0) s)))
		 |t|)))
	 (if (zerop prob-val) ; 確率0の場合は最小のdouble-floatの値を返す
	     least-positive-double-float
	     prob-val)))
     mu Sigma Sigma-inv W s <<1>>-acc Lambda~-acc)))

;;; posteriori distribution of NGnet
;; posteriori requires list of joint distribution
(defun posteriori-distribution (joint-lst)
  (if (= (length joint-lst) 1)
      '(1.0d0) ; 1 unitの場合
      (let ((denominator (ssum joint-lst)))
	(mapcar (lambda (x)
		  (let ((prob-val (/ x denominator)))
		    (if (zerop prob-val) ; 確率0の場合は最小のdouble-floatの値を返す
			least-positive-double-float
			prob-val)))
		joint-lst))))

;;; weighted average (online version)
;; accumulator は統計量(スカラー,ベクトル,行列)のリスト
(defun weighted-average (x y func accumulator eta posteriori)
  (mapcar
   (lambda (acc post)
     (if (arrayp acc)
	 (m+ (m* (- 1 eta) acc)
	     (m* eta (funcall func x y) post))
	 (+ (* (- 1 eta) acc)
	    (* eta (funcall func x y) post))))
   accumulator posteriori))

;;; Accumulators

(defun <<1>> (x y acc eta posteriori)
  (weighted-average x y
		    (lambda (x y)
		      (declare (ignore x y))
		      1)
		    acc eta posteriori))

(defun <<x>> (x y acc eta posteriori)
  (weighted-average x y
		    (lambda (x y)
		      (declare (ignore y))
		      x)
		    acc eta posteriori))

(defun <<x*x^T>> (x y acc eta posteriori)
  (weighted-average x y
		    (lambda (x y)
		      (declare (ignore y))
		      (if (arrayp x)
			  (m* x (m-t x))
			  (* x x)))
		    acc eta posteriori))

;; for output standard deviation

(defun <<y^T*y>> (x y acc eta posteriori)
  (weighted-average x y
		    (lambda (x y)
		      (declare (ignore x))
		      (if (arrayp y)
			  (aref (m* (m-t y) y) 0 0)
			  (* y y)))
		    acc eta posteriori))

(defun <<x~*y^T>> (x y acc eta posteriori)
  (weighted-average x y
		    (lambda (x y)
		      (let ((x~ (ext-vector x)))
			(if (arrayp y)
			    (m* x~ (m-t y))
			    (m* y x~))))
		    acc eta posteriori))

;; for output liner regression unit W

(defun <<y*x~^T>> (x y acc eta posteriori)
  (weighted-average x y
		    (lambda (x y)
		      (let ((x~ (ext-vector x)))
			(if (arrayp y)
			    (m* y x~)
			    (m* y (m-t x~))
			    )))
		    acc eta posteriori))

(defun <<x~*x~^T>> (x y acc eta posteriori)
  (weighted-average x y
		    (lambda (x y)
		      (declare (ignore y))
		      (let ((x~ (ext-vector x)))
			(m* x~ (m-t x~))))
		    acc eta posteriori))

;; (4.6)
(defun Lambda~ (x acc eta posteriori)
  (mapcar
   (lambda (a.Lambda~ post)
     (let ((x~ (ext-vector x)))
       (m* (/ 1 (- 1 eta))
	   (m- a.Lambda~
	       (m* (/ 1d0 (+ (/ 1d0 eta) -1d0
			     (* post (aref (m* (m-t x~) a.Lambda~ x~) 0 0))))
		   post a.Lambda~ x~ (m-t x~) a.Lambda~)))))
   acc posteriori))

;; (4.5a)
(defun calc-mu (<<1>>-acc <<x>>-acc)
  (mapcar (lambda (a.1 a.x)
	    (if (arrayp a.x)
		(m* (/ 1.0d0 a.1) a.x)		
		(/ a.x a.1)))
	  <<1>>-acc <<x>>-acc))

;; (4.5b)
(defun calc-Sigma (<<1>>-acc <<x*x^T>>-acc mu)
  (mapcar (lambda (a.1 a.x*x^T a.mu)
	    (if (arrayp a.x*x^T)
		(m- (m* (/ 1.0d0 a.1) a.x*x^T)
		    (m* a.mu (m-t a.mu)))
		(- (/ a.x*x^T a.1) (* a.mu a.mu))))
	  <<1>>-acc <<x*x^T>>-acc mu))

;; (5.6) regularized on-line update of Sigma
(defun regularization-of-Sigma (Sigma alpha)
  (let ((N (array-dimension Sigma 0)))
    (m+ Sigma (m* (* alpha (/ (trace-matrix Sigma) N)) (umat N)))))

(defun calc-Sigma-reg (<<1>>-acc <<x*x^T>>-acc <<delta^2>>-acc mu alpha)
  (let ((N (array-dimension (car mu) 0)))
    (mapcar (lambda (a.1 a.x*x^T a.delta^2 a.mu)
	      (m* (/ 1.0d0 a.1)
		  (m+ a.x*x^T
		      (m* (- a.1) a.mu (m-t a.mu))
		      (m* (* alpha a.delta^2) (umat N)))))
	    <<1>>-acc <<x*x^T>>-acc <<delta^2>>-acc mu)))

;; (4.7) N*N partial matrix of Lambda~ correspond to inverse matrix of Sigma
(defun calc-Sigma-invert (<<1>>-acc Lambda~-acc)
  (let ((N (1- (array-dimension (car Lambda~-acc) 0))))
    (mapcar (lambda (a.1 a.Lambda~)	    
	      (partial-matrix (m* a.Lambda~ a.1) N N))
	    <<1>>-acc Lambda~-acc)))

;; (4.8b)			      
(defun calc-W (x y posteriori Lambda~-acc eta W)
  (let ((x~ (ext-vector x)))
    (mapcar (lambda (post a.Lambda~ a.W)
	      (m+ a.W
		  (m* eta post
		      (if (arrayp y)
			  (m- y (m* a.W x~))
			  (- y (aref (m* a.W x~) 0 0)))
		      (m-t x~)
		      a.Lambda~)))
	    posteriori Lambda~-acc W)))

;; (4.5d)
(defun calc-s (<<1>>-acc <<y^T*y>>-acc <<x~*y^T>>-acc W)
  (let ((D (array-dimension (car W) 0)))
    (mapcar (lambda (a.1 a.y^T*y a.x~*y^T a.W)
	      (abs (/ (- a.y^T*y (trace-matrix (m* a.W a.x~*y^T)))
		      (* a.1 D))))
	    <<1>>-acc <<y^T*y>>-acc <<x~*y^T>>-acc W)))
;; (4.2b)
(defun update-eta (eta-pre lambda-factor)
  (/ 1 (+ 1 (/ lambda-factor eta-pre))))

;;; 5.2 Regularization of on-line EM algorithm

;; (5.10)
(defun calc-eta-post-delta^2 (x mu.old mu.new eta post <<1>>-acc)
  (let ((N (array-dimension (car mu.new) 0)))
    (mapcar (lambda (a.mu.old a.mu.new a.post a.1)
	      (let ((x-mu.new (m- x a.mu.new))
		    (mu.new-mu.old (m- a.mu.new a.mu.old)))
		(/ (+ (* eta a.post (square (euclidean-norm x-mu.new)))
		      (* (- 1 eta) (square (euclidean-norm mu.new-mu.old)) a.1))
		   N)))
	    mu.old mu.new post <<1>>-acc)))

;; (5.9)
(defun <<delta^2>> (eta-post-delta^2 acc eta)
  (mapcar (lambda (a.eta-post-delta^2 a.delta^2)
	    (+ a.delta^2 a.eta-post-delta^2 (- (* eta a.delta^2))))
	  eta-post-delta^2 acc))

;; (5.12)
(defun regularize-Lambda~ (Lambda~ eta-post-delta^2 eta posteriori alpha)
  (let ((N (1- (array-dimension (car Lambda~) 0))))
    (mapcar (lambda (a.Lambda~ a.eta-post-delta^2 a.post)
	      (nlet itr ((L a.Lambda~)
			 (i 0))
		(if (= i N)
		    L
		    (let ((nu-coefficient (sqrt (* alpha (/ a.eta-post-delta^2 eta a.post))))
			  (nu (make-vector (1+ N) :initial-element 0.0d0)))
		      (setf (aref nu i 0) nu-coefficient)
		      (itr (m- L (m* (/ 1.0d0 (1+ (aref (m* (* eta a.post) (m-t nu) L nu) 0 0)))
				    (* eta a.post) L nu (m-t nu) L))
			   (1+ i))))))
	    Lambda~ eta-post-delta^2 posteriori)))

;;; 再帰処理のときに受け渡すデータを構造体sessionとして定義
;;; sessionはパラメータと，<<>>で囲われている再帰計算の途中経過で渡される統計量を保持する
(defstruct session
  mu Sigma Sigma-inv W s
  <<1>>-acc              ;dim: 1
  <<x>>-acc           ;dim: N
  <<x*x^T>>-acc    ;dim: N*N
  <<y*x~^T>>-acc   ;dim: D*(N+1)
  <<x~*x~^T>>-acc    ;dim: (N+1)*(N+1)
  Lambda~-acc            ;dim: (N+1)*(N+1)
  <<y^T*y>>-acc  ;dim: 1
  <<x~*y^T>>-acc   ;dim: (N+1)*D
  <<delta^2>>-acc        ;dim: 1
  (eta 0.5d0) (lambda-factor 0.99d0)
  (alpha 0.1d0) ; 正則化パラメータ
  (LL 0d0) ; 対数尤度 nilだと対数尤度の記録、表示は行なわない
  (global-count 1))

(defun session-params (session)
  (list (session-mu session) (session-Sigma session) (session-Sigma-inv session) (session-W session) (session-s session) (session-<<1>>-acc session) (session-Lambda~-acc session)))

(defun oem-1step (x y session)
  (let* ((eta (session-eta session))
	 (joint (apply #'joint-probability x y (session-params session)))
	 (post (posteriori-distribution joint))
	 (<<1>>-acc (<<1>> x y (session-<<1>>-acc session) eta post))
	 (<<x>>-acc (<<x>> x y (session-<<x>>-acc session) eta post))
	 (mu (calc-mu <<1>>-acc <<x>>-acc)) ; mu update
	 (eta-post-delta^2-lst (if (arrayp x) (calc-eta-post-delta^2 x (session-mu session)
								     mu eta post <<1>>-acc)))
	 (<<delta^2>>-acc (if (arrayp x) (<<delta^2>> eta-post-delta^2-lst
						      (session-<<delta^2>>-acc session) eta)))
	 (<<x*x^T>>-acc (<<x*x^T>> x y (session-<<x*x^T>>-acc session) eta post))
	 (Sigma (if (arrayp x)
		    (calc-Sigma-reg <<1>>-acc <<x*x^T>>-acc <<delta^2>>-acc mu (session-alpha session))
		    (calc-Sigma <<1>>-acc <<x*x^T>>-acc mu))) ; Sigma update
	 ;(Sigma-ori (calc-Sigma <<1>>-acc <<x*x^T>>-acc mu))
	 (<<y*x~^T>>-acc (<<y*x~^T>> x y (session-<<y*x~^t>>-acc session) eta post))
	 (<<x~*x~^T>>-acc (<<x~*x~^T>> x y (session-<<x~*x~^t>>-acc session) eta post))
	 (Lambda~-acc (let ((Lambda~-acc (Lambda~ x (session-Lambda~-acc session) eta post)))
			(if (arrayp x)
			    (regularize-Lambda~ Lambda~-acc eta-post-delta^2-lst eta post (session-alpha session))
			    Lambda~-acc)))
	 (Sigma-inv (calc-Sigma-invert <<1>>-acc Lambda~-acc))
	 (W (calc-W x y post Lambda~-acc eta (session-W session))) ; W update
	 (<<y^T*y>>-acc (<<y^T*y>> x y (session-<<y^T*y>>-acc session) eta post))
	 (<<x~*y^T>>-acc (<<x~*y^T>> x y (session-<<x~*y^T>>-acc session) eta post))
	 (s (calc-s <<1>>-acc <<y^T*y>>-acc <<x~*y^T>>-acc W)) ; s update
	 )
    ;; (plot-prediction-1dim *last-session* (seq (- (* 1.5 pi)) (* 1.5 pi) :by 0.05))
    ;; (plot-joint-1dim *last-session* ;learned-session
    ;; 		     (seq (- (* 1.5 pi)) (* 1.5 pi) :by 0.1)
    ;; 		     (seq (- 2.0) 2.0 :by 0.1))
    
    ;; output log-likelihood (if LL slot is null, it will not output this value)
    (if (session-LL session) (format t "LL= ~A~%" (/ (session-LL session) (session-global-count session))))
    (finish-output)
    (make-session :mu              mu
		  :Sigma           Sigma
		  :Sigma-inv       Sigma-inv
		  :W               W
		  :s               s
		  :<<1>>-acc       <<1>>-acc
		  :<<x>>-acc       <<x>>-acc
		  :<<x*x^T>>-acc   <<x*x^T>>-acc
		  :<<y*x~^T>>-acc  <<y*x~^T>>-acc
		  :<<x~*x~^T>>-acc <<x~*x~^T>>-acc
		  :Lambda~-acc     Lambda~-acc
		  :<<y^T*y>>-acc   <<y^T*y>>-acc
		  :<<x~*y^T>>-acc  <<x~*y^T>>-acc
		  :<<delta^2>>-acc <<delta^2>>-acc
		  :eta             (update-eta eta (session-lambda-factor session))
		  :lambda-factor   (session-lambda-factor session)
		  :alpha           (session-alpha session)
		  :LL              (if (session-LL session) (+ (log (ssum joint)) (session-LL session)))
		  :global-count    (1+ (session-global-count session)))))

(defvar *last-session*)

;;; シミュレーションのメインループ
;; データがリストで与えられているとき用
(defun oem-data-list (x.dat y.dat session)
  (if (null x.dat)
      session
      (let ((sub-session (oem-1step (car x.dat) (car y.dat) session)))
	(setf *last-session* sub-session)
	(oem-data-list (cdr x.dat) (cdr y.dat) sub-session))))

;;; evaluation for dataset
;;; log-likelihood of dataset 対数尤度関数
(defun log-likelihood (x.dat y.dat parameters)
  (ssum
   (mapcar
    (lambda (x y)
      (let ((joint-sum (ssum (apply #'joint-probability x y parameters))))
	(if (zerop joint-sum)
	    0d0
	    (log joint-sum))))
    x.dat y.dat)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; 予測関連 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; 予測線
(defun prediction (x session)
  (let* ((x~ (ext-vector x))
	 (G (mapcar (lambda (mu-i Sigma-i Sigma-inv-i)
		      (universal-gaussian x mu-i Sigma-i)) ; :inv-Sigma Sigma-inv-i
		    (session-mu session) (session-Sigma session) (session-Sigma-inv session)))
	 (sum-G (ssum G))
	 (N (mapcar (lambda (Gi) (/ Gi sum-G)) G))
	 (W (session-W session)))
    (if (= (array-dimension (car W) 0) 1)
	(ssum (mapcar (lambda (Ni Wi) (aref (m* Ni Wi x~) 0 0)) N W))
	(ssum-vec (mapcar (lambda (Ni Wi) (m* Ni Wi x~)) N W)))))

;;; 事後分布のリストからユニットIDを選択
(defun unit-selector (post)
  (let ((prob-sum-list
	 (nlet itr ((post-lst post)
		     (sum 0)
		     (product '()))
	   (if (= (length post-lst) 1)
	       (reverse (cons 1.0 product))
	       (let ((prob-sum (+ sum (car post-lst))))
		 (itr (cdr post-lst) prob-sum (cons prob-sum product)))))))    
    (position-if (lambda (x)
		   (< (random 1.0d0) x)) prob-sum-list)))

;;; 出力の分布p(y | x, t, theta)からのサンプル
(defun sample-from-output-distribution (x unit-id session)
  (let* ((Wi (nth unit-id (session-W session)))
	 (si (nth unit-id (session-s session)))
	 (x-vec (ext-vector x))
	 (mean-vec (m* Wi x-vec)))
    (if (= (array-dimension Wi 0) 1)	
	(random-normal :mean (aref mean-vec 0 0) :sd si)
	(list->vector (loop for d from 0 to (1- (array-dimension Wi 0))
			 collect (random-normal :mean (aref mean-vec d 0) :sd si))))))

(defun predict-sampling (x session)
  (let* ((joint (apply #'joint-probability
		       x
		       (prediction x session)
		       (session-params session)))
	 (post (posteriori-distribution joint))
	 (unit-id (unit-selector post)))
    (sample-from-output-distribution x unit-id session)))

;;; データセットが与えられているとき
;; あるデータセットを頭から走査し,お尻まで達したら再び頭からアクセスする
;; x.datとy.datは最初に設定されて以降変わらないデータそのもの．
;; stateはバッファであり,これをnullになるまでcdrしていく
;; 一つのデータセットに対して何度も繰り返してオンラインに学習する際に使う
(defclass data-stream ()
  ((x.dat    :accessor x.dat-of    :initform () :initarg :x.dat)
   (y.dat   :accessor y.dat-of   :initform () :initarg :y.dat)
   (x.state  :accessor x.state-of  :initform () :initarg :x.state)
   (y.state :accessor y.state-of :initform () :initarg :y.state)))

(defun make-data-stream (x.dat y.dat)
  (let ((d (make-instance 'data-stream :x.dat x.dat :y.dat y.dat)))
    (setf (x.state-of d) x.dat)
    (setf (y.state-of d) y.dat)
    d))

;; return: xとyの1個ずつのデータを多値で返す
(defmethod read-from-data-stream ((d data-stream))
  (if (null (x.state-of d)) 
      (progn ; stateが空のとき
	(setf (x.state-of  d) (cdr (x.dat-of  d)))
	(setf (y.state-of d) (cdr (y.dat-of d)))
	(values (car (x.dat-of d)) (car (y.dat-of d))))
      (let ((x  (car (x.state-of  d)))
	    (y (car (y.state-of d))))
	(setf (x.state-of  d) (cdr (x.state-of  d)))
	(setf (y.state-of d) (cdr (y.state-of d)))
	(values x y))))

(defmethod clear-data-stream! ((d data-stream))
  (let ((d d))
    (setf (x.state-of  d) (cdr (x.dat-of  d)))
    (setf (y.state-of d) (cdr (y.dat-of d)))
    d))

;;; データセットを生成する関数が与えられているとき
;;; スロットに入出力の対応を与える関数を持たせる
(defclass data-stream-general-source ()  
  (;; 入力を生成する関数. 引数なし. 通常ある区間での一様乱数になる.
   (input-generator :accessor input-generator-of
		    :initarg :input-generator)
   ;; 入力を引数にとって出力を返す関数. 入力の形式はinput-generatorに従う
   (input-output-function :accessor input-output-function-of
			  :initarg :input-output-function)))

(defmethod read-from-data-stream ((d data-stream-general-source))
  (let* ((x  (funcall (input-generator-of d)))
	 (y (funcall (input-output-function-of d) x)))
    (values x y)))

;; データストリームで与えられているとき用
;; データリストを循環して学習し続ける
(defun oem-data-stream (d-stream session)
  (multiple-value-bind (x y)
      (read-from-data-stream d-stream)
    (let ((sub-session (oem-1step x y session)))
      (setf *last-session* sub-session)
      (oem-data-stream d-stream sub-session))))

;; 一様分布から2次元配列の要素をサンプリングしてくる
(defun make-random-array (dim-list low-val high-val)
  (let ((a (make-array dim-list)))
    (sfor (i 0 (1- (car dim-list)))
      (sfor (j 0 (1- (cadr dim-list)))
	(setf (aref a i j) (random-uniform low-val high-val))))
    a))

;;; パラメータを初期化したsession構造体を返す
;; 一様乱数の範囲などは適宜調整が必要. アキュームレータの初期値は大体このままでもよい.
;; 収束判定の期間は計算時間や最終的な精度に関して重要になる
;; LLをnilにすると対数尤度の出力をやめる
(defun make-initialized-session (n-of-unit input-dimension output-dimension
			    &key (lambda-factor 0.99d0) (alpha 0.1d0) (LL 0d0))
  (make-session
   :mu (n-times-collect n-of-unit
	 (if (= input-dimension 1)
	     (random-uniform -1.0d0 1.0d0)
	     (make-random-array (list input-dimension 1) -1.0d0 1.0d0)))
   :Sigma (n-times-collect n-of-unit
	    (if (= input-dimension 1)
		1d0
		(umat input-dimension)))
   :Sigma-inv (n-times-collect n-of-unit
		(if (= input-dimension 1)
		    1d0
		    (umat input-dimension)))
   :W (n-times-collect n-of-unit
	(make-random-array (list output-dimension (1+ input-dimension)) -1.0d0 1.0d0))
      
   :s (n-times-collect n-of-unit (random-uniform 0.1d0 1.0d0))

   :<<1>>-acc (make-list n-of-unit :initial-element 0.1d0)
   :<<x>>-acc  (make-list n-of-unit :initial-element ; dim: N
			     (if (= input-dimension 1)
				 0.1d0
				 (make-vector input-dimension :initial-element 0.1d0)))
   :<<x*x^T>>-acc (make-list n-of-unit :initial-element ; dim: N*N
				   (if (= input-dimension 1)
				       0.1d0
				       (make-array (list input-dimension input-dimension)
						   :initial-element 0.1d0)))
   ;; :<<x-mu>>-acc (make-list n-of-unit :initial-element ; dim: N*N
   ;; 			       (if (= input-dimension 1)
   ;; 				   0.1d0
   ;; 				   (make-array (list input-dimension input-dimension)
   ;; 					       :initial-element 0.1d0)))
   
   :<<y*x~^T>>-acc (make-list n-of-unit :initial-element ; dim: D*(N+1)
				    (make-array (list output-dimension (1+ input-dimension))
						:initial-element 0.1d0))
   :<<x~*x~^T>>-acc (make-list n-of-unit ; dim: (N+1)*(N+1)
				   :initial-element (umat (1+ input-dimension)))
   :Lambda~-acc (make-list n-of-unit ; dim: (N+1)*(N+1)
   			   :initial-element (umat (1+ input-dimension)))
   
   :<<y^T*y>>-acc (make-list n-of-unit :initial-element 0.1d0)
   :<<x~*y^T>>-acc (make-list n-of-unit :initial-element ; dim: (N+1)*D
				    (make-array (list (1+ input-dimension) output-dimension)
						:initial-element 0.1d0))
   :<<delta^2>>-acc (make-list n-of-unit :initial-element 0.1d0)
   :eta 0.5d0
   :lambda-factor lambda-factor
   :alpha alpha
   :LL LL))
