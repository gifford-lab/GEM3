
(define (makeBBG bayes prob size)
  (BayesBindingGenerator. bayes
			  prob
			  size
			  Boolean.FALSE$))
(define (makeBPG bayes prob size probmin sizemin)
  (BayesBindingGenerator. bayes
			  prob
			  size
			  probmin
			  sizemin))


(define (makemethod name bindinggenerator)
  (list name bindinggenerator))
(define methodname car)
(define methodbg cadr)

(define (genenamesFromFile inputport)
  (let ((factor (read inputport))
	(target (read inputport)))
    (if (or (eof-object? factor)
	    (eof-object? target))
	'()
	(cons target (genenamesFromFile inputport)))))

(define (genesFromFile inputport)
  (let ((genenames (genenamesFromFile inputport)))
    (filter (lambda (x) (not (null? x)))
	    (map (lambda (x) (tryCatch( .getGene genome (symbol->string x))
				      (lambda (e) ;(.printStackTrace e) 
					'())))
		 genenames))))

(define (regionsToFile iter outputport)
  (if (.hasNext iter)
      (begin (display (.toString (.next iter))
		    outputport)
	     (newline outputport)
	     (regionsToFile iter outputport))))

(define (regionsListToFile lst outputport)
  (map (lambda (x) (begin (display (.toString x) outputport)
			  (newline outputport)))
       lst))

(define (vizregions regions)
  (visualize vizargs 
	     (MapperIterator. (RegionExpander. 2000 2000)
			      (if (list? regions)
				  (set! regions (lstToIterator regions))
				  regions))))

(define (reportDA distances)
  (if (> (length distances) 0)
      (let* ((nomatch (length (filter (lambda (x) (= x -1)) distances)))
	     (toofar (length (filter (lambda (x) (= x -2)) distances)))
	     (n (length distances))
	     (med (median (filter (lambda (x) (>= x 0)) distances)))
	     (avg (mean (filter (lambda (x) (>= x 0)) distances)))
	     (variance (mean (map (lambda (x) (expt (- x avg) 2)) distances))))
	(string-append ": mean=" avg "  median="  med  "  n=" n "  variance=" variance "  meanvariance=" (/ variance n) "  meanstddev=" (sqrt (/ variance n)) "  nomatch=" nomatch "  toofar=" toofar))
      " no distances "))

(define (distanceAnalysis source sink regionsize distlimit regioniterator )
  (let* ((regions (MapperIterator. (RegionExpander. (/ regionsize 2)
						    (/ regionsize 2))
				   (ExpanderIterator. source regioniterator)))
	 (distances (ExpanderIterator. (DistanceGenerator. source sink)
				       regioniterator)))
    distances))

(define (doDA regions motifs method)  
  (let* ((distances (iteratorTolst (distanceAnalysis (methodbg method) motifs 1000 5000 regions))))
    (display (string-append (methodname method)
			    (reportDA distances)
			    "\n"))
    distances))

; motif centric version that comes up with one average distance per motif (
;  average across all binding events nearest to that motif) and then reports
; statistics on those averages
(define (doMotifDA regions motifs method)
  (let* ((mtbb (MotifToBinding. motifs (methodbg method)))
	 (distances (iteratorTolst (ExpanderIterator. mtbb regions))))
    (string-append (methodname method)
		   (reportDA distances)
		   "\n")))
	 
    


;;
;; does curve that gives % binding included by distance d from motif
;;

(define (motifBindingCoverage motifs method regiter)
  (let* ((mtbb (MotifToBinding. motifs (methodbg method))))
    (.coverage mtbb regiter)))

(define (reportCoverage coverage)
  (let* ((distances '(10 20 30 40 50 60 70 80 90 100))
	 (total (Array.get coverage  (- (Array.getLength coverage) 1)))
	 (dists (map (lambda (d) (Array.get coverage d)) distances)))
    (string-append (map (lambda (x) (trunc (/ x total) 2))
			dists)
		   " " (Array.get coverage (- (Array.getLength coverage) 2))
		   " " (Array.get coverage (- (Array.getLength coverage) 1))
		   )))


;;
;;  useful tools for testing thresholds to pick the best one to use
;;

(define (filterbe events t1 t2)
  (filter (lambda (event) (and (>= (.getConf event) t1)
			       (>= (.getSize event) t2)))
	  events))

(define (bindingInRegions method regions)
  (let ((regionList (iteratorTolst regions)))
    (iteratorTolst (ExpanderIterator. (methodbg method)
				      (lstToIterator regionList)))))

(define (bindingInFileRegions method filename)
  (bindingInRegions  method (MapperIterator. (Gene2IGRRegion. genome)
					     (lstToIterator (genesFromFile (open-input-file filename))))))

(define maxfp 2)
(define maxwg 10000)
(define (test-thresholds-versions loosestmethod thresholds accept?)  
  (let* ((fpresults (bindingInRegions loosestmethod (negregions)))
	 (wgresults (bindingInRegions loosestmethod (wholegenome)))
	 (okthresholds (filter (lambda (thresh)
				 (and 
				  (<= (length (filter (lambda (be) (accept? thresh be)) fpresults))
				      maxfp)
				  (<= (length (filter (lambda (be) (accept? thresh be)) wgresults))
				      maxwg)))
				 
			       thresholds)))
;    (display (accumulate (lambda (x y)
;			   (string-append x "\n" y)) 
;			 ""
;			 (map (lambda (thresh)
;				(let ((n (length (filter (lambda (be) (accept? thresh be)) fpresults))))
;				  (list thresh n)))
;			      thresholds)))
    (letrec ((findbest (lambda (maxfound maxthresh rest)
			 (if (null? rest)
			     (list maxfound maxthresh)
			     (let ((thisfound (length (filter (lambda (be) (accept? (car rest)
										    be))
							      wgresults))))
;			       (display (string-append "thresh " (car rest)
;						       " numfound " thisfound
;						       "\n"))
			       (if (>= thisfound maxfound)
				   (findbest thisfound (car rest) (cdr rest))
				   (findbest maxfound maxthresh (cdr rest))))))))
      (findbest -1 '(none none) okthresholds))))

;;
;; ROC 
;;

(define (countBoundRegions bg regionsiter)
  (let* ((regions (iteratorTolst regionsiter))
	 (numbound (accumulate + 0 (map (lambda (r)
					  (let ((iter (.execute bg r)))
					    (if (.hasNext iter)
						1
						0)))
					regions))))
    (list numbound (length regions))))

(define (getSensSpec bg)  
  (let ((bg (if (list? bg) (methodbg bg) bg)))
    (let ((posinfo (countBoundRegions bg (posregions)))
	  (neginfo (countBoundRegions bg (negregions))))
      (let ((sens (/ (* 1.0 (car posinfo)) (cadr posinfo)))
	    (spec (- 1 (/ (* 1.0 (car neginfo)) (cadr neginfo)))))
	(list sens spec posinfo neginfo )))))


(define (makeVisualizer args)
  (import "edu.mit.csail.psrg.viz.*")
  (import "edu.mit.csail.psrg.viz.sauron.*")
  (let* ((opts (GFFSVOptions.))
	 (opts (GFFStreamVisualizer.parseArgs (lstToArray args) opts))
	 (visualizer (GFFStreamVisualizer. opts)))
    visualizer))

(define verylargedistance 10000000000)
(define (closest region others)
  (define (h others best)
    (if (null? others)
	best
	(let ((d (.distance region (car others))))
	  (h (cdr others)
	     (if (< d best) d best)))))
  (h others verylargedistance))