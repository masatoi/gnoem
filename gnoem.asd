(defsystem "gnoem"
    :description "GNOEM: Normalized Gaussian Network with On-line EM algorithm."
    :version "0.01"
    :author "Satoshi Imai <satoshi.imai@gmail.com>"
    :licence "Public Domain"
    :depends-on (:wiz-util)
    :components ((:file "packages")
		 (:file "matrix")
		 (:file "multivariate-gaussian" :depends-on ("matrix"))
		 (:file "gnoem" :depends-on ("multivariate-gaussian"))))
