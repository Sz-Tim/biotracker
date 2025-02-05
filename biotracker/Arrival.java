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
    private final String sourceSiteID;
    private final String arrivalSiteID;
    private final int arrivalStatus;
    private final double arrivalDepth;
    private final double arrivalDensity;

    public Arrival(String sourceSiteID, String arrivalSiteID, int arrivalStatus, double arrivalDepth, double arrivalDensity) {
        this.sourceSiteID = sourceSiteID;
        this.arrivalSiteID = arrivalSiteID;
        this.arrivalStatus = arrivalStatus;
        this.arrivalDepth = arrivalDepth;
        this.arrivalDensity = arrivalDensity;
    }

    public String getSourceSiteID() {
        return this.sourceSiteID;
    }

    public String getArrivalSiteID() {
        return this.arrivalSiteID;
    }

    public int getArrivalStatus() {
        return this.arrivalStatus;
    }

    public double getArrivalDepth() {
        return this.arrivalDepth;
    }

    public double getArrivalDensity() {
        return this.arrivalDensity;
    }
}
