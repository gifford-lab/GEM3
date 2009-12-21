/*
 * Author: tdanford
 * Date: Nov 10, 2008
 */
package edu.mit.csail.cgs.utils;

public interface BitVector {

	/**
	 * Gives a hex representation of the BitVector.  Leading bits (i.e., if the length of the
	 * vector is not divisible by 4) are assumed to be zeros.
	 * @return The String object containing the Hex representation.
	 */
	public String toHexString();

	public boolean isOn(int index);
	public boolean isOff(int index);

	public int length();

	public int countOnBits();

	public void turnOnBit(int index);
	public void turnOffBit(int index);
	
	public static class Complete implements BitVector { 
		private int length;
		private boolean value;
		
		public Complete(int len, boolean v) { 
			length = len; 
			value = v;
		}

		public int countOnBits() {
			return value ? length : 0;
		}

		public boolean isOff(int index) {
			return !value;
		}

		public boolean isOn(int index) {
			return value;
		}

		public int length() {
			return length;
		}

		public String toHexString() {
			int fs = length/4;
			StringBuilder sb = new StringBuilder();
			if(value) {
				for(int i = 0; i < fs; i++) { 
					sb.append("F");
				}
				int extra = length%4;
				int value = 1;
				for(int i = 0, base = 8; i < extra; i++, base /= 2) { 
					value += base; 
				}
				
				if(value < 10) { 
					sb.append(String.valueOf(value)); 
				} else { 
					int c = (char)('A' + (value-10));
					sb.append(c);
				}
				
			} else { 
				for(int i = 0; i <= fs; i++) { 
					sb.append("0");
				}
			}
			return sb.toString();
		}

		public void turnOffBit(int index) {
			throw new UnsupportedOperationException(String.format("BitVector.Complete.turnOffBit()"));
		}

		public void turnOnBit(int index) {
			throw new UnsupportedOperationException(String.format("BitVector.Complete.turnOnBit()"));
		}
	}
}