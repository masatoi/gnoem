;;; -*- Coding: utf-8; Mode: Lisp -*-

(defpackage :wiz-numerical
  (:use :common-lisp :wiz-util)
  (:export :matrixp :num-rows :num-cols :square-matrix? :make-matrix
	   :make-identity-matrix :copy-matrix :print-matrix
	   :transpose-matrix :multiply-matrix :add-matrix :subtract-matrix
	   :invert-matrix :solve-matrix :trace-matrix :partial-matrix
	   ;; arrayed-vector and matrix wrapper function
	   :make-vector :list->vector :simple-vector->arrayed-vector
	   :vector-cat :vecot-cat2 :vector-length
	   :euclidean-norm :m* :m+ :ssum-vec :m- :m-t :umat :zero-mat :m-1
	   :m-append-horizon :vec :mat :mapmat :mapvec :cholesky :det-cholesky
	   ;; gaussian
	   :gaussian :multivariate-gaussian :universal-gaussian :random-uniform :random-normal))

(defpackage :gnoem
  (:use :common-lisp :wiz-util :wiz-numerical)
  (:export :oem-1step :oem-data-list :oem-data-stream
	   :prediction :make-data-stream :read-from-data-stream
	   :make-initialized-session))
