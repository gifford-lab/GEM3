/*
 * Author: tdanford
 * Date: Dec 3, 2008
 */
package edu.mit.csail.cgs.utils.models.bns;

import java.util.*;

import edu.mit.csail.cgs.utils.models.Model;
import edu.mit.csail.cgs.utils.models.data.DataFrame;

public class BNTest {
	
	public static void main(String[] args) { 
		BNTest test = new BNTest(1000);
		BN<ObsModel> bn = new BN<ObsModel>(test.obs, "a", "b", "c", "d", "e");
		
		BNSearch searcher = new BNSearch(bn);
		searcher.search(100);

		System.out.println("Final Network:");
		bn.print();
		
		for(int i = 0; i < 10; i++) { 
			System.out.println(bn.sample().toString());
		}
		
		ObsModel testObs = test.sampleObs();
		System.out.println("\nTest: " + testObs.toString());
		System.out.println("b Posterior: " + bn.posterior(testObs, "b"));
	}
	
	public DataFrame<ObsModel> obs;
	private Random rand;
	
	public BNTest(int n) {
		rand = new Random();
		obs = new DataFrame<ObsModel>(ObsModel.class);
		
		for(int i = 0; i < n; i++) { 
			obs.addObject(sampleObs());
			//System.out.println(m.toString());
		}
	}
	
	public ObsModel sampleObs() { 
		int fiveten = rand.nextInt(2);
		String c = fiveten == 0 ? "five" : "ten";
		Integer a = rand.nextInt(3);
		Integer b = (fiveten == 0 ? a + 5 : a + 10) + (rand.nextInt(3)-1);
		Integer d = rand.nextInt(5);
		String e  = d % 2 == 0 ? "even" : "odd";
		ObsModel m = new ObsModel(a, b, c, d, e);
		return m;
	}

	public static class ObsModel extends Model { 
		public Integer a, b;
		public String c;
		public Integer d; 
		public String e;
		
		public ObsModel() { 
			super();
		}
		
		public ObsModel(Integer a, Integer b, String c, Integer d, String e) { 
			this.a = a;
			this.b = b;
			this.c = c;
			this.d = d;
			this.e = e;
		}
	}
}
