/*
 * Created on Apr 12, 2007
 */
package edu.mit.csail.cgs.echo.gui;

import edu.mit.csail.cgs.ewok.types.SelfDescribingConstant;
import edu.mit.csail.cgs.utils.Factory;

public interface SelfDescribingConstantFactory extends Factory<SelfDescribingConstant> {
    public boolean isImmediate();
}
