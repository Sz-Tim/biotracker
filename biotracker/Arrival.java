/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotracker;

/**
 * This class is used to record an instance of an arrival within the particle class,
 * so that particles can build a list of arrival instances that they make, independently
 * of the simulation-wide arrays. This avoid concurrent updating, but requires postprocessing
 * at the end of the run to count arrivals.
 * <p>
 * NOT PRESENTLY USED 21/11/18
 *
 * @author sa01ta
 */
public class Arrival {
    private final int sourceLocation;
    private final int arrivalLocation;
    private final double arrivalTime;
    private final double arrivalDensity;

    public Arrival(int sourceLocation, int arrivalLocation, double arrivalTime, double arrivalDensity) {
        this.sourceLocation = sourceLocation;
        this.arrivalLocation = arrivalLocation;
        this.arrivalTime = arrivalTime;
        this.arrivalDensity = arrivalDensity;
    }

    public int getSourceLocation() {
        return this.sourceLocation;
    }

    public int getArrivalLocation() {
        return this.arrivalLocation;
    }

    public double getArrivalTime() {
        return this.arrivalTime;
    }

    public double getArrivalDensity() {
        return this.arrivalDensity;
    }
}
