;;; -*- coding: utf-8 mode: lisp -*-

;; GNOEMをロード
(asdf:oos 'asdf:load-op :gnoem)

(in-package :gnoem)

;; 近似対象の関数
(defun schaal-function (x1 x2)
  (max (exp (* -10d0 x1 x1))
       (exp (* -50d0 x2 x2))
       (* 1.25d0 (exp (* -5d0 (+ (* x1 x1)
				 (* x2 x2)))))))

;; 入出力データを1000個用意
(defparameter x-data
  (loop repeat 1000 collect
       (make-vector 2 :initial-contents (list (1- (random 2d0)) (1- (random 2d0))))))

(defparameter y-data
  (mapcar (lambda (x_in)
	    (schaal-function (aref x_in 0 0) (aref x_in 1 0)))
	  x-data))

;; 新規にセッションをつくる
;; (make-initialized-session ユニット数 入力次元 出力次元)
(defparameter init-session
  (make-initialized-session 300 2 1 :lambda-factor 0.99d0 :alpha 0.1d0))

;; リストに格納されているデータから学習
(defparameter learned-session (oem-data-list x-data y-data init-session))

;; 予測
(prediction (make-vector 2 :initial-contents (list 0.5d0 0.5d0)) learned-session)
;; => 0.14885842827023427d0
