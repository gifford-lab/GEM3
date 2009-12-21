; imports
(import "edu.mit.csail.psrg.ewok.*")
(import "edu.mit.csail.psrg.ewok.nouns.*")
(import "edu.mit.csail.psrg.ewok.verbs.*")
(import "java.util.*")
(import "javax.swing.*")
(import "java.io.*")
(import "edu.mit.csail.psrg.Datasets.species.*")

; general setup to make life pleasant
(define nil '())
; java List to scheme list
(define (filter pred lst)
  (define (helper left sofar)
    (if (null? left)
	sofar
	(if (pred (car left))
	    (helper (cdr left)
		    (cons (car left)
			  sofar))
	    (helper (cdr left)
		    sofar))))
  (helper lst '()))

(define (accumulate proc init lst)
  (if (null? lst)
      init
      (accumulate proc (proc (car lst) init) (cdr lst))))

(define (firstn lst n)
  (if (or (null? lst)
	  (<= n 0))
      '()
      (cons (car lst)
	    (firstn (cdr lst)
		    (- n 1)))))

(define (trunc number digits)
  (let ((ex (expt 10 digits)))
    (/ (round (* ex number)) ex)))

(define-method (getList (lst List))
  (let ((output nil)
	(lislen (.size lst)))    
    (define (r pos)
      (if (= pos lislen)
	  '()
	  (cons (.get lst pos)
		(r (+ pos 1)))))
    (r 0)))
; java Iterator to scheme list
(define-method (iteratorTolst (iter Iterator))
  (define (helper sofar)
    (if (.hasNext iter)	
	(helper (append sofar (list (.next iter))))
	sofar))
  (helper '()))

;scheme list to java Iterator
(define (lstToIterator lst)
  (.iterator (makeList lst)))
;scheme list to java List
(define (makeList lst)
  (let ((javalist (ArrayList.)))
    (letrec ((helper (lambda (x)
		    (if (null? x)
			javalist
			(begin (.add javalist (car x))
			       (helper (cdr x)))))))
      (helper lst))))
(define (ListTolst lst)
  (letrec ((helper (lambda (i output)
		     (if (= i (.size lst))
			 output
			 (helper (+ i 1)
				 (append output (list (.get lst i))))))))
    (helper 0 '())))

(define (lstToArray lst)
  (letrec ((array (Array.newInstance String.class (length lst)))
	   (helper (lambda (i left)
		     (if (null? left)
			 array
			 (begin 
			   (Array.set array i (car left))
			   (helper (+ i 1) (cdr left)))))))
    (helper 0 lst)))

(define (ArrayTolst array) 
  (letrec ((helper (lambda (i output)
		     (if (= i (Array.getLength array))
			 output
			 (helper (+ i 1)
				 (append output (list (Array.get array i))))))))
    (helper 0 '())))
		   
(define (lstToSet lst)
  (let ((javaset (HashSet.)))
    (letrec ((helper (lambda (x)
		       (if (null? x)
			   javaset
			   (begin (.add javaset (car x))
				  (helper (cdr x)))))))
      (helper lst))))

(define (uniq lst)
  (iteratorTolst (.iterator (lstToSet lst))))
(define (uniqiter . iters)
  (define (recurse iter set)
    (if (.hasNext iter)
	(begin 
	  (.add set (.next iter))
	  (recurse iter set))
	set))
  (let ((set (HashSet.)))
    (map (lambda (i) 
	   (recurse i set))
	 iters)
    (.iterator set)))

(define (mean lst)
  (define (recurse sum count left)
    (if (null? left)
	(if (= count 0)
	    0
	    (/ sum count))
	(recurse (+ sum (car left))
		 (+ count 1)
		 (cdr left))))
  (recurse 0 0 lst))

(define (median lstOfNums)  
  (let ((javaList (ArrayList.)))
    (for-each (lambda(x) (.add javaList (Double. (.toString x)))) lstOfNums)
    (let ((array (.toArray javaList)))
      (Arrays.sort array)
      (Array.get array (/ (Array.getLength array) 2)))))

(define (variance lst)
  (let ((xbar (mean lst)))
    (mean (map (lambda (x) (expt (- x xbar) 2)) lst))))


(define (commaseparate lst)
  (accumulate (lambda (a b) (if b
				(if a
				    (string-append a ", " b)
				    b)
				a))
	      #f
	      lst))
