package biotracker;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.filefilter.WildcardFileFilter;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.nio.file.FileSystems;
import java.sql.Array;
import java.util.*;
import java.util.concurrent.*;
import java.util.stream.IntStream;

/**
 * @author tomdude
 */
public class Particle_track {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws Exception {

        System.out.println("Starting biotracker\n");
        Date date = new Date();
        System.out.println(date);

        long heapMaxSize = Runtime.getRuntime().maxMemory();
        System.out.println("Max heap " + heapMaxSize);

        long startTime = System.currentTimeMillis();

        System.out.println("Reading in data\n");
        System.out.println(System.getProperty("user.dir"));

        RunProperties rp = new RunProperties(args[0]); // first (and only?) cmd line arg is properties filename e.g. model_setup.properties
        rp.checkForConflictingProperties();

        ISO_datestr currentIsoDate = new ISO_datestr(rp.start_ymd);
        ISO_datestr endIsoDate = new ISO_datestr(rp.end_ymd);

        int numberOfDays = endIsoDate.getDateNum() - currentIsoDate.getDateNum() + 1;

        // Print all main arguments
        System.out.println(rp.printProperties());

        // --------------------------------------------------------------------------------------
        // File reading and domain configuration
        // --------------------------------------------------------------------------------------       
        List<Mesh> meshes = new ArrayList<>(2);
        meshes.add(new Mesh(rp.mesh0, rp.meshType0, rp.coordOS));
        if (!rp.mesh1.isEmpty()) {
            meshes.add(new Mesh(rp.mesh1, rp.meshType1, rp.coordOS));
            if (rp.hfFilePrefix1.equalsIgnoreCase("westcoms2")) {
                int adjElem1 = switch (rp.hfFilePrefix0) {
                    case "etive28" -> 69784;
                    default -> 69784;
                };
                int adjElem0 = switch (rp.hfFilePrefix0) {
                    case "etive28" -> 16;
                    default -> 16;
                };
                meshes.get(0).setAdjoiningElement(adjElem0);
                meshes.get(1).setAdjoiningElement(adjElem1);
            }
        }
        int[] allelems = IntStream.rangeClosed(0, meshes.get(0).getUvnode()[0].length - 1).toArray();

        double subStepDt = rp.dt / (double) rp.stepsPerStep; // number of seconds per substep
        double dev_perstep = Math.pow(0.1, subStepDt);
        System.out.println("Particle subStepDt = " + subStepDt + " dev_perstep = " + dev_perstep);

        // --------------------------------------------------------------------------------------
        // Creating initial particle array
        // --------------------------------------------------------------------------------------
        List<HabitatSite> habitat;
        List<HabitatSite> habitatEnd;
        System.out.println("Creating start and end sites");
        habitat = IOUtils.createHabitatSites(rp.sitefile, null, 4, false, meshes, rp);
        habitatEnd = IOUtils.createHabitatSites(rp.sitefileEnd, null, 4, false, meshes, rp);

        List<String> siteStartNames = new ArrayList<>();
        for (HabitatSite habitatSite : habitat) {
            siteStartNames.add(habitatSite.getID());
        }
        List<String> siteEndNames = new ArrayList<>();
        for (HabitatSite habitatSite : habitatEnd) {
            siteEndNames.add(habitatSite.getID());
        }

        // --------------------------------------------------------------------------------------
        // Starting density of each particle
        // --------------------------------------------------------------------------------------
        double startDensity = 1.0;
        double[][] siteDensities = new double[habitat.size()][numberOfDays];
        if(!rp.siteDensityPath.isEmpty()) {
            siteDensities = IOUtils.readDailyDensities(rp.siteDensityPath, siteStartNames, rp);
        } else {
            for (int i = 0; i < habitat.size(); i++) {
                for (int j = 0; j < numberOfDays; j++) {
                    siteDensities[i][j] = startDensity;
                }
            }
        }

        // --------------------------------------------------------------------------------------
        // Setup particles
        // --------------------------------------------------------------------------------------
        int nparts = rp.nparts * habitat.size();
        List<Particle> particles = new ArrayList<>(nparts);
        int numParticlesCreated = 0; // Counter to keep track of how many particles have been created
        boolean allowRelease = true; // boolean to be switched after a single release event

        // --------------------------------------------------------------------------------------
        // Read in particles to be restarted, if there are any
        // --------------------------------------------------------------------------------------
        if (!rp.restartParticles.equalsIgnoreCase("")) {
            List<Particle> restartParts = IOUtils.readRestartParticles(rp);
            particles.addAll(restartParts);
            numParticlesCreated = numParticlesCreated + (restartParts.size());
            System.out.println("numberOfParticles: " + numParticlesCreated + " " + particles.size());
        }

        // --------------------------------------------------------------------------------------
        // Setup hydrodynamic fields and file lists
        // --------------------------------------------------------------------------------------
        List<HydroField> hydroFields = new ArrayList<>(Collections.nCopies(3, null));
        List<HydroField> hfToday = new ArrayList<>(Collections.nCopies(3, null));
        List<HydroField> hfTomorrow = new ArrayList<>(Collections.nCopies(3, null));
        List<Integer> hfIters = new ArrayList<>();
        ArrayList<ArrayList<String>> missingHydroFiles = new ArrayList<>(Collections.nCopies(3, null));
        for (int m = 0; m < meshes.size(); m++) {
            hfIters.add(m);
            missingHydroFiles.set(m, IOUtils.checkHydroFilesExist(rp, currentIsoDate, endIsoDate, numberOfDays, m));
            if (!missingHydroFiles.get(m).isEmpty()) {
                System.out.println("\nWarning! Cannot find these files for mesh " + m + ":");
                for (String missing: missingHydroFiles.get(m)) {
                    System.out.println(missing);
                }
                System.out.println("-- The previous day will be used for each instead!");
            } else {
                System.out.println("  All found for mesh " + m);
            }
        }
        if (rp.needStokes) {
            int m = 2;
            hfIters.add(m);
            missingHydroFiles.set(m, IOUtils.checkHydroFilesExist(rp, currentIsoDate, endIsoDate, numberOfDays, m));
            if (!missingHydroFiles.get(m).isEmpty()) {
                System.out.println("\nWarning! Cannot find these files for mesh " + m + ":");
                for (String missing: missingHydroFiles.get(m)) {
                    System.out.println(missing);
                }
                System.out.println("-- The previous day will be used for each instead!");
            } else {
                System.out.println("  All found for mesh " + m);
            }
        }

        // --------------------------------------------------------------------------------------
        // Load daylight hours
        // --------------------------------------------------------------------------------------
        int[][] daylightHours = new int[numberOfDays][2];
        boolean isDaytime = false;
        if (!rp.daylightPath.isEmpty()) {
            daylightHours = IOUtils.readDaylightHours(rp.daylightPath, currentIsoDate, endIsoDate, numberOfDays, " ", false);
        }

        // --------------------------------------------------------------------------------------
        // Initialize counters
        // --------------------------------------------------------------------------------------
        int stepcount = 0;
        double elapsedHours = 0; // updated in HOURS as the simulation progresses
        int[] freeViableSettleExit;
        int[] particleStatus = new int[]{0,0,0,0,0,0};

        // --------------------------------------------------------------------------------------
        // Parallel details
        // --------------------------------------------------------------------------------------
        int numberOfExecutorThreads = rp.parallelThreads;
        System.out.println("Number of executor threads = " + numberOfExecutorThreads);
        ExecutorService executorService = Executors.newFixedThreadPool(numberOfExecutorThreads);
        CompletionService<List<Particle>> executorCompletionService = new ExecutorCompletionService<>(executorService);

        // --------------------------------------------------------------------------------------
        // Initialize storage objects
        // --------------------------------------------------------------------------------------
        final Collection<Callable<List<Particle>>> callables = new ArrayList<>();

        // Set up arrays to hold particle density*hour counts
        int pstepsInd2 = 1;
        if (rp.splitPsteps) {
            pstepsInd2 = habitat.size();
        }
        float[][] pstepsImmature = new float[meshes.get(0).getNElems()][pstepsInd2];
        float[][] pstepsMature = new float[meshes.get(0).getNElems()][pstepsInd2];

        // Set up array to hold connectivity counts
        float[][] connectivity = new float[habitat.size()][habitatEnd.size()];
        float[][] connectivityDepth2 = new float[habitat.size()][habitatEnd.size()];
        float[][] connectivityImmature = new float[habitat.size()][habitatEnd.size()];
        float[][] connectivityImmatureDepth2 = new float[habitat.size()][habitatEnd.size()];

        int[][] elemActivity = new int[meshes.get(0).getNElems()][3]; // count of sink, swim, float within each element
        int[][] hourActivity = new int[numberOfDays*24+1][3]; // count of sink, swim, float within each hour
        float[][] vertDistrImmature = new float[meshes.get(0).getNElems()][rp.vertDistrMax+1]; // bin upper z limits: from 0 to vertDistrMax
        float[][] vertDistrMature = new float[meshes.get(0).getNElems()][rp.vertDistrMax+1]; // bin upper z limits: from 0 to vertDistrMax

        String particleRestartHeader = "hour,ID,startDate,age,startLocation,x,y,elem,status,density,mesh,depth,depthLayer,degreeDays,xTot,yTot,xyTot,zTot";
        String arrivalHeader = "ID,startDate,startTime,startLocation,endDate,endTime,endLocation,age,density";

        if (rp.recordMovement) {
            IOUtils.writeMovementsHeader("ID,date,hour,step,startDate,age,density,x,y,z,layer,mesh,status,degreeDays,sink,swim,temp,salinity,mortality,tempSurface,dX,dY,dZ",
                    "movementFile.csv");
        }

        try {
            // --------------------------------------------------------------------------------------
            // Start time loop
            // --------------------------------------------------------------------------------------
            long currTime = System.currentTimeMillis();
            for (int fnum = 0; fnum < numberOfDays; fnum++) {

                // ******* Is it the last day of the simulation? *******
                // If so, readHydroFields will read the last day file twice, and use the first hour of endDay as hour1
                // of endDay+1 for interpolation purposes.
                // Other design choices possible:
                // - stop particles at hour 23 on last day (introduces and hour gap if using for restart)
                // - do no interpolation for last hour of run
                // Whatever you choose - a small error would be introduced. Could do using actual next day file but this
                // means losing a day's worth of hydrodynamic data if running in operational mode.
                boolean isLastDay = (fnum == numberOfDays - 1) && rp.duplicateLastDay;

                String today = currentIsoDate.getDateStr();
                if (rp.recordLocations) {
                    IOUtils.printFileHeader(particleRestartHeader, "locations_" + today + ".csv");
                }
                if (rp.recordArrivals) {
                    IOUtils.printFileHeader(arrivalHeader, "arrivals_" + today + ".csv");
                }

                long splitTime = System.currentTimeMillis();
                System.out.printf("----------- Day %d (%s) - Stepcount %d (%.1f hrs) ----------- \n",
                        fnum + 1, currentIsoDate, stepcount, elapsedHours);
                System.out.printf("Elapsed time:     %.1f h (%.0f s)\n", (splitTime - startTime) / 1000.0 / 3600, (splitTime - startTime) / 1000.0);
                System.out.printf("Last 24hr time:   %.1f s\n", (splitTime - currTime) / 1000.0);
                System.out.println("-------------------");
                currTime = System.currentTimeMillis();

                // default, run loop forwards
                // ---- LOOP OVER ENTRIES IN THE HYDRO OUTPUT ------------------------
                for (int currentHour = 0; currentHour < 24; currentHour++) {

                    // Calculate current time of the day (complete hours elapsed since midnight)
                    if (!rp.daylightPath.isEmpty()) {
                        isDaytime = daylightHours[fnum][0] <= currentHour && daylightHours[fnum][1] >= currentHour;
                    }

                    // Read new hydrodynamic fields?
                    if (currentHour == 0) {
                        String currentDateShort = currentIsoDate.toStringShort();
                        ISO_datestr nextIsoDate = ISO_datestr.getTomorrow(currentIsoDate);
                        String nextDateShort = nextIsoDate.toStringShort();
                        for (Integer h : hfIters) {
                            int m = h==2 ? 0 : h;
                            boolean missingToday = missingHydroFiles.get(h).contains(currentDateShort);
                            boolean missingTomorrow = missingHydroFiles.get(h).contains(nextDateShort);
                            if (missingToday) {
                                System.out.print("Can't find file: Reusing hydrodynamics from yesterday\n");
                                System.out.println("-------------------");
                            } else {
                                hfToday.set(h, fnum == 0 ? readHydroField(meshes.get(m), currentIsoDate, rp, h) : hfTomorrow.get(h));
                                hfTomorrow.set(h, (missingTomorrow | isLastDay) ? hfToday.get(h) : readHydroField(meshes.get(m), nextIsoDate, rp, h));
                                if (h==2) {
                                    hydroFields.set(h, mergeHydroFieldsStokes(hfToday.get(h), hfTomorrow.get(h), rp));
                                } else {
                                    hydroFields.set(h, mergeHydroFields(hfToday.get(h), hfTomorrow.get(h), rp));
                                }
                                System.out.println("-------------------");
                            }
                        }
                    }

                    if (currentHour == 12) {
                        System.out.println();
                    }
                    System.out.printf("%02d", currentHour);

                    for (int i=0; i < habitat.size(); i++) {
                        int siteElem = habitat.get(i).getContainingFVCOMElem();
                        int m = habitat.get(i).getContainingMesh();
                        int nLayers = (int) Mesh.findNearestSigmas(20.0, meshes.get(m).getSiglay(), (float) habitat.get(i).getDepth())[0][0];
                        double[][] currentConditions = new double[nLayers][13];
                        double[] siteLoc = new double[2];
                        siteLoc[0] = habitat.get(i).getLocation()[0];
                        siteLoc[1] = habitat.get(i).getLocation()[1];
                        double[] siteNodeDist = Particle.distanceToNodes(siteLoc, siteElem, meshes.get(0).getNodexy(), meshes.get(0).getTrinodes(), rp.coordOS);
                        double[] siteStokesUV = {0,0};
                        for (int j = 0; j < nLayers; j++) {
                            currentConditions[j][0] = hydroFields.get(m).getU()[currentHour][j][siteElem];
                            currentConditions[j][1] = hydroFields.get(m).getV()[currentHour][j][siteElem];
                            currentConditions[j][2] = hydroFields.get(m).getW()[currentHour][j][siteElem];
                            currentConditions[j][3] = Math.sqrt(currentConditions[j][0]*currentConditions[j][0] + currentConditions[j][1]*currentConditions[j][1]);
                            currentConditions[j][4] = rp.needS ? hydroFields.get(m).getAvgFromTrinodes(meshes.get(m), siteLoc, j, siteElem, currentHour, "salinity", rp) : -9999;
                            currentConditions[j][5] = rp.needT ? hydroFields.get(m).getAvgFromTrinodes(meshes.get(m), siteLoc, j, siteElem, currentHour, "temp", rp) : -9999;
                            currentConditions[j][6] = rp.needK ? hydroFields.get(m).getAvgFromTrinodes(meshes.get(m), siteLoc, j, siteElem, currentHour, "km", rp) : -9999;
                            if (rp.needStokes) {
                                siteStokesUV = hydroFields.get(2).getAvgFromTrinodes(siteElem, meshes.get(0).getTrinodes(), siteNodeDist, currentHour, 0);
                            }
                            currentConditions[j][7] = rp.needStokes ? siteStokesUV[0] : -9999;
                            currentConditions[j][8] = rp.needStokes ? siteStokesUV[1] : -9999;
                            currentConditions[j][9] = rp.needStokes ? Math.sqrt(siteStokesUV[0]*siteStokesUV[0] + siteStokesUV[1]*siteStokesUV[1]) : -9999;
                            currentConditions[j][10] = rp.needStokes ? hydroFields.get(2).getAvgFromTrinodes(meshes.get(m), siteLoc, 0, siteElem, currentHour, "Hsig", rp) : -9999;
                            currentConditions[j][11] = rp.needStokes ? hydroFields.get(2).getAvgFromTrinodes(meshes.get(m), siteLoc, 0, siteElem, currentHour, "Dir", rp) : -9999;
                            currentConditions[j][12] = rp.needStokes ? hydroFields.get(2).getAvgFromTrinodes(meshes.get(m), siteLoc, 0, siteElem, currentHour, "Tm", rp) : -9999;
                            habitat.get(i).addEnvCondition(currentConditions[j]);
                        }
                        // calculate eggs per female per timestep based on temperature and scale initial particle densities
                        // These are nearly identical under the temperatures in Scotland (r = 0.998)
                        float[][] eggSigmas = Mesh.findNearestSigmas(10, meshes.get(m).getSiglay(), (float) habitat.get(i).getDepth());
                        double eggTemperature = hydroFields.get(m).getValueAtDepth(meshes.get(m), siteElem, siteLoc, 2, currentHour, "temp", rp, eggSigmas);
                        double eggsPerFemale = switch (rp.eggTemp_fn) {
                            case "constant" -> rp.eggTemp_b.get(0) * (rp.dt / (60 * 60 * 24));
                            case "linear" ->
                                    (rp.eggTemp_b.get(0) + rp.eggTemp_b.get(1) * eggTemperature) * (rp.dt / (60 * 60 * 24));
                            case "quadratic" ->
                                    (rp.eggTemp_b.get(0) * Math.pow((eggTemperature + rp.eggTemp_b.get(1)), 2)) * (rp.dt / (60 * 60 * 24));
                            case "logistic" ->
                                    (rp.eggTemp_b.get(0) / (1 + Math.exp(-rp.eggTemp_b.get(1) * (eggTemperature - rp.eggTemp_b.get(2)))) + rp.eggTemp_b.get(3)) * (rp.dt / (60 * 60 * 24));
                            default -> rp.eggTemp_b.get(0) * (rp.dt / (60 * 60 * 24));
                        };
                        habitat.get(i).setScale(siteDensities[i][fnum] / rp.nparts * eggsPerFemale);
                    }

                    // Create new particles, if releases are scheduled hourly, or if release is scheduled for this hour
                    double hourRemainder = Math.round((elapsedHours % rp.releaseInterval) * 100.0) / 100.0;
                    if ((rp.releaseScenario == 0 && elapsedHours >= rp.releaseTime && allowRelease) ||
                            rp.releaseScenario == 1 && (hourRemainder == 0 || hourRemainder == rp.releaseInterval) ||
                            (rp.releaseScenario == 2 && elapsedHours >= rp.releaseTime && elapsedHours <= rp.releaseTimeEnd)) {
                        List<Particle> newParts = createNewParticles(habitat, meshes, rp, currentIsoDate, currentHour, numParticlesCreated);
                        particles.addAll(newParts);
                        numParticlesCreated += newParts.size();
                        // If only one release to be made, prevent further releases
                        if (rp.releaseScenario == 0) {
                            allowRelease = false;
                        }
                        System.out.print("+");
                    } else {
                        System.out.print(".");
                    }

                    // ---- INTERPOLATE BETWEEN ENTRIES IN THE HYDRO OUTPUT ------------------------
                    for (int step = 0; step < rp.stepsPerStep; step++) {

                        // MOVE the particles
                        if (rp.parallelThreads > 1) {
                            int particlesSize = particles.size();
                            int listStep = particlesSize / numberOfExecutorThreads;
                            for (int i = 0; i < numberOfExecutorThreads; i++) {
                                List<Particle> subList;
                                if (i == numberOfExecutorThreads - 1) {
                                    // Note: ArrayList.subList(a,b) is inclusive of a, but exclusive of b =>
                                    subList = particles.subList(i * listStep, particlesSize);
                                } else {
                                    subList = particles.subList(i * listStep, (i + 1) * listStep);
                                }
                                callables.add(new ParallelParticleMover(subList, elapsedHours, currentHour, step, subStepDt, rp,
                                        meshes, hydroFields, habitatEnd, allelems, isDaytime, today, elemActivity, hourActivity));
                            }
                            for (Callable<List<Particle>> callable : callables) {
                                executorCompletionService.submit(callable);
                            }
                            for (Callable<List<Particle>> callable : callables) {
                                executorCompletionService.take().get();
                            }
                            callables.clear();
                        } else {
                            for (Particle part : particles) {
                                ParallelParticleMover.move(part, elapsedHours, currentHour, step, subStepDt, rp,
                                        meshes, hydroFields, habitatEnd, allelems, isDaytime, today, elemActivity, hourActivity);
                            }
                        }

                        elapsedHours += subStepDt / 3600.0;
                        stepcount++;
                    }

                    // COUNT the number of particles in different states
                    particleStatus = particleCountStatus(particles, particleStatus);

                    if (rp.recordLocations) {
                        IOUtils.particlesToRestartFile(particles, currentHour, "locations_" + today + ".csv", true, rp, 1); // rp.nparts * rp.numberOfDays * 10
                    }

                    // update connectivity
                    if (rp.recordConnectivity) {
                        for (Particle particle : particles) {
                            if (particle.getArrivals() != null) {
                                for (Arrival arrival : particle.getArrivals()) {
                                    int destIndex = siteEndNames.indexOf(arrival.getArrivalSiteID());
                                    if (arrival.getArrivalStatus() == 2) { // copepodids
                                        if (arrival.getArrivalDepth() >= rp.connectDepth1_min &&
                                                arrival.getArrivalDepth() <= rp.connectDepth1_max) {
                                            connectivity[particle.getStartIndex()][destIndex] += (float) arrival.getArrivalDensity();
                                        }
                                        if (rp.recordConnectivityDepth2 &&
                                                arrival.getArrivalDepth() >= rp.connectDepth2_min &&
                                                arrival.getArrivalDepth() <= rp.connectDepth2_max) {
                                            connectivityDepth2[particle.getStartIndex()][destIndex] += (float) arrival.getArrivalDensity();
                                        }
                                    } else if (rp.connectImmature && arrival.getArrivalStatus() == 1) {
                                        if (arrival.getArrivalDepth() >= rp.connectDepth1_min &&
                                                arrival.getArrivalDepth() <= rp.connectDepth1_max) {
                                            connectivityImmature[particle.getStartIndex()][destIndex] += (float) arrival.getArrivalDensity();
                                        }
                                        if (rp.recordConnectivityDepth2 &&
                                                arrival.getArrivalDepth() >= rp.connectDepth2_min &&
                                                arrival.getArrivalDepth() <= rp.connectDepth2_max) {
                                            connectivityImmatureDepth2[particle.getStartIndex()][destIndex] += (float) arrival.getArrivalDensity();
                                        }
                                    }
                                }
                            }
                        }
                    }

                    // Hourly updates to pstep arrays
                    if (rp.recordPsteps) {
                        IOUtils.pstepsUpdater(particles, rp, pstepsMature, pstepsImmature, rp.dt);
                    }

                    // Write Psteps
                    if (rp.recordPsteps && stepcount % (rp.pstepsInterval * rp.stepsPerStep) == 0) {
                        if (rp.recordImmature) {
                            IOUtils.writeNonZeros2DArrayToCSV(pstepsImmature, "i,NA,value","%d,%d,%.4e" , "pstepsImmature_" + today + "_" + Math.round(elapsedHours) + ".csv");
                            pstepsImmature = new float[meshes.get(0).getNElems()][pstepsInd2];
                        }
                        IOUtils.writeNonZeros2DArrayToCSV(pstepsMature, "i,NA,value","%d,%d,%.4e", "pstepsMature_" + today + "_" + Math.round(elapsedHours) + ".csv");
                        pstepsMature = new float[meshes.get(0).getNElems()][pstepsInd2];
                    }

                    if (rp.recordConnectivity && stepcount % (rp.connectivityInterval * rp.stepsPerStep) == 0) {
                        IOUtils.writeNonZeros2DArrayToCSV(connectivity, "source,destination,value", "%d,%d,%.4e",
                                "connectivity_" + rp.connectDepth1_min + "-" + rp.connectDepth1_max + "m_" + today + "_" + Math.round(elapsedHours) + ".csv");
                        connectivity = new float[habitat.size()][habitatEnd.size()];
                        if (rp.recordConnectivityDepth2) {
                            IOUtils.writeNonZeros2DArrayToCSV(connectivityDepth2, "source,destination,value", "%d,%d,%.4e",
                                    "connectivity_" + rp.connectDepth2_min + "-" + rp.connectDepth2_max + "m_" + today + "_" + Math.round(elapsedHours) + ".csv");
                            connectivityDepth2 = new float[habitat.size()][habitatEnd.size()];
                        }
                        if (rp.connectImmature) {
                            IOUtils.writeNonZeros2DArrayToCSV(connectivityImmature, "source,destination,value", "%d,%d,%.4e",
                                    "connectivityImm_" + rp.connectDepth1_min + "-" + rp.connectDepth1_max + "m_" + today + "_" + Math.round(elapsedHours) + ".csv");
                            connectivityImmature = new float[habitat.size()][habitatEnd.size()];
                            if (rp.recordConnectivityDepth2) {
                                IOUtils.writeNonZeros2DArrayToCSV(connectivityImmatureDepth2, "source,destination,value", "%d,%d,%.4e",
                                        "connectivityImm_" + rp.connectDepth2_min + "-" + rp.connectDepth2_max + "m_" + today + "_" + Math.round(elapsedHours) + ".csv");
                                connectivityImmatureDepth2 = new float[habitat.size()][habitatEnd.size()];
                            }
                        }
                        for (Particle particle : particles) {
                            particle.clearArrivals();
                        }
                    }

                    if (rp.recordVertDistr) {
                        IOUtils.vertDistrUpdater(particles, rp, vertDistrMature, vertDistrImmature, rp.dt);
                        if (stepcount % (rp.vertDistrInterval * rp.stepsPerStep) == 0) {
                            if (rp.recordImmature) {
                                IOUtils.writeNonZeros2DArrayToCSV(vertDistrImmature, "i,z,value", "%d,%d,%.4e", "vertDistrImmature_" + today + "_" + Math.round(elapsedHours) + ".csv");
                                vertDistrImmature = new float[meshes.get(0).getNElems()][rp.vertDistrMax+1];
                            }
                            IOUtils.writeNonZeros2DArrayToCSV(vertDistrMature, "i,z,value", "%d,%d,%.4e", "vertDistrMature_" + today + "_" + Math.round(elapsedHours) + ".csv");
                            vertDistrMature = new float[meshes.get(0).getNElems()][rp.vertDistrMax+1];
                        }
                    }

                    // Clean up "dead" (666) and "exited" (66) particles
                    particles.removeIf(part -> part.getStatus() == 666 || part.getStatus() == 66);
                }
                System.out.println();
                System.out.println("-------------------");
                System.out.printf("Total particles:  %,d\n", numParticlesCreated);
                System.out.printf("Nauplii:          %,d\n", particleStatus[1]);
                System.out.printf("Copepodids:       %,d\n", particleStatus[2]);
                System.out.printf("Exited domain:    %,d\n", particleStatus[4]);
                System.out.printf("Dead:             %,d\n", particleStatus[5]);
                System.out.println("-------------------");

                FileWriter fstream = new FileWriter("siteConditions_" + today + "_" + Math.round(elapsedHours) + ".csv", false);
                PrintWriter out = new PrintWriter(fstream);
                out.println("site,x,y,depth,mesh,centroid,elem,MeshType," +
                        "u_Avg,v_Avg,w_Avg,uv_Avg,salinity_Avg,temperature_Avg,k_Avg,stokesU0_Avg,stokesV0_Avg,stokesUV0_Avg,Hsig_Avg,waveDir_Avg,waveT_Avg," +
                        "u,v,w,uv,salinity,temperature,k,stokesU0,stokesV0,stokesUV0,Hsig,waveDir,waveT");
                for (HabitatSite habitatSite : habitat) {
                    out.println(habitatSite);
                    habitatSite.setEnvConditionDay(new double[13]);
                    habitatSite.setEnvConditionCountDay(new int[13]);
                }
                out.close();

                System.out.println();

                if (rp.backwards) {
                    currentIsoDate.takeDay();
                } else {
                    currentIsoDate.addDay();
                }
            }

            // Write out the final locations of the particles.
            // Note that the last hour of the last day has by now been iterated over, and the day has been advanced
            // to the day after the simulation finished.
            // So this is the location of the particles at t=0 on the day after the last simulated day, ready to 
            // start a new run on the next day.
            IOUtils.printFileHeader(particleRestartHeader, "locationsEnd_" + currentIsoDate.getDateStr() + ".csv");
            IOUtils.particlesToRestartFile(particles, 0, "locationsEnd_" + currentIsoDate.getDateStr() + ".csv", true, rp, 1);
            if (rp.recordActivity) {
                IOUtils.writeIntegerArrayToFile(elemActivity, "elementActivity.csv");
                IOUtils.writeIntegerArrayToFile(hourActivity, "hourActivity.csv");
            }

            executorService.shutdownNow();
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            executorService.shutdownNow();
        }

        FileWriter fstream = new FileWriter("startSitesUsed.csv", false);
        PrintWriter out = new PrintWriter(fstream);
        out.println("site,x,y,depth,mesh,centroid,elem,MeshType," +
                "u_Avg,v_Avg,w_Avg,uv_Avg,salinity_Avg,temperature_Avg,k_Avg,stokesU0_Avg,stokesV0_Avg,stokesUV0_Avg,Hsig_Avg,waveDir_Avg,waveT_Avg," +
                "u,v,w,uv,salinity,temperature,k,stokesU0,stokesV0,stokesUV0,Hsig,waveDir,waveT");
        for (HabitatSite habitatSite : habitat) {
            out.println(habitatSite);
        }
        out.close();

        long endTime = System.currentTimeMillis();
        System.out.println("Elapsed time = " + (endTime - startTime) / 1000.0);
    }

    /**
     * Method to create new particles. These must be appended to the existing list
     *
     * @return List of the new particles to be appended to existing list
     */
    public static List<Particle> createNewParticles(List<HabitatSite> habitat, List<Mesh> meshes, RunProperties rp,
                                                    ISO_datestr currentDate, int currentTime, int numParticlesCreated) {

        List<Particle> newParts = new ArrayList<>(rp.nparts * habitat.size());
        for (int i = 0; i < rp.nparts * habitat.size(); i++) {
            int startid = i % habitat.size();
            double xstart = habitat.get(startid).getLocation()[0];
            double ystart = habitat.get(startid).getLocation()[1];
            int meshStart = habitat.get(startid).getContainingMesh();
            int elemFVCOMStart = habitat.get(startid).getContainingFVCOMElem();
            int[] elemROMSStartU = habitat.get(startid).getContainingROMSElemU();
            int[] elemROMSStartV = habitat.get(startid).getContainingROMSElemV();
            int[] nearestROMSGridPointU = habitat.get(startid).getNearestROMSPointU();
            int[] nearestROMSGridPointV = habitat.get(startid).getNearestROMSPointV();
            double startDensity = habitat.get(startid).getScale();

            Particle p = new Particle(xstart, ystart, rp.startDepth, habitat.get(startid).getID(), startid, numParticlesCreated + i,
                    startDensity, currentDate, currentTime, rp.coordOS);
            p.setMesh(meshStart);
            p.setElem(elemFVCOMStart);
            p.setLayerFromDepth(meshes.get(meshStart).getDepthUvnode()[elemFVCOMStart], meshes.get(meshStart).getSiglay());
            p.setROMSElemU(elemROMSStartU);
            p.setROMSElemV(elemROMSStartV);
            p.setROMSnearestPointU(nearestROMSGridPointU);
            p.setROMSnearestPointV(nearestROMSGridPointV);

            newParts.add(p);
        }
        return newParts;
    }

    /**
     * Method to handle the various cases of reading in hydrodynamic model output files.
     * The cases handled are:
     * i) Single mesh, FVCOM only. In this case, the current day is read in in its entirety,
     * plus the first hour of the next day.
     * ii) More than one mesh, hour 0-22 of the day. In this case, two hours of data are read for
     * each relevant mesh. For an FVCOM mesh, this is just two hours from the same .nc file.
     * For a ROMS mesh, this is two separate files.
     * iii) More than one mesh, hour 23 of the day. In this case, as single hour (23) is read from
     * the first file, and then record 0 from tomorrow's file is read.
     */
//    public static List<HydroField> readHydroFields(List<Mesh> meshes, ISO_datestr currentIsoDate, int currentHour, boolean isLastDay, RunProperties rp) throws InterruptedException {
//        List<HydroField> hydroFields = new ArrayList<>();
//
//        // 24 hr files only case - read once a day
//        int m = 0;
//        for (Mesh mesh : meshes) {
//            if (mesh.getType().equalsIgnoreCase("FVCOM")) {
//
//                // Occasionally the read fails for unknown reasons when the .nc file exists, so try up to 5 times before quitting
//                int readAttempt = 0;
//                while (readAttempt < 5) {
//                    try{
//                        // Dima file naming format: minch2_20171229_0003.nc
//                        String[] varNames1 = new String[]{"u", "v", "ww", "salinity", "temp", "zeta", "kh", "viscofh", "short_wave", "Hsig", "Dir", "Tm01"};
//
//                        // Normal "forwards time"
//                        if (!rp.backwards) {
//                            ISO_datestr tomorrow = ISO_datestr.getTomorrow(currentIsoDate);
//                            if (isLastDay) {
//                                System.out.println("** Last day or missing day - reading same hydro file twice **");
//                                tomorrow = currentIsoDate;
//                            }
//                            // Inelegant, but it works and I'm in a hurry. Clean it up if you're procrastinating on something else.
//                            String[] dirsT1 = new String[]{
//                                    rp.datadir + rp.datadirPrefix + currentIsoDate.getYear() + rp.datadirSuffix + FileSystems.getDefault().getSeparator(),
//                                    rp.datadir2 + rp.datadir2Prefix + currentIsoDate.getYear() + rp.datadir2Suffix + FileSystems.getDefault().getSeparator()};
//                            String[] dirsT2 = new String[]{
//                                    rp.datadir + rp.datadirPrefix + tomorrow.getYear() + rp.datadirSuffix + FileSystems.getDefault().getSeparator(),
//                                    rp.datadir2 + rp.datadir2Prefix + tomorrow.getYear() + rp.datadir2Suffix + FileSystems.getDefault().getSeparator()};
//                            String[] filesT1 = new String[]{
//                                    rp.mesh1Domain + "_" + currentIsoDate.getYear() + String.format("%02d", currentIsoDate.getMonth()) + String.format("%02d", currentIsoDate.getDay()) + "*.nc",
//                                    rp.mesh2Domain + "_" + currentIsoDate.getYear() + String.format("%02d", currentIsoDate.getMonth()) + String.format("%02d", currentIsoDate.getDay()) + "*.nc"};
//                            String[] filesT2 = new String[]{
//                                    rp.mesh1Domain + "_" + tomorrow.getYear() + String.format("%02d", tomorrow.getMonth()) + String.format("%02d", tomorrow.getDay()) + "*.nc",
//                                    rp.mesh2Domain + "_" + tomorrow.getYear() + String.format("%02d", tomorrow.getMonth()) + String.format("%02d", tomorrow.getDay()) + "*.nc"};
//                            List<File> files1 = (List<File>) FileUtils.listFiles(new File(dirsT1[m]), new WildcardFileFilter(filesT1[m]), null);
//                            List<File> files2 = (List<File>) FileUtils.listFiles(new File(dirsT2[m]), new WildcardFileFilter(filesT2[m]), null);
//                            // Read both files and combine
//                            hydroFields.add(new HydroField(files1.get(0).getCanonicalPath(), files2.get(0).getCanonicalPath(), varNames1, null, null, null, "FVCOM", rp));
//                        }
//                        // Instead read time backwards, so need yesterday instead
//                        else {
//                            List<File> files1 = (List<File>) FileUtils.listFiles(
//                                    new File(rp.datadir + rp.datadirPrefix + currentIsoDate.getYear() + rp.datadirSuffix + FileSystems.getDefault().getSeparator()),
//                                    new WildcardFileFilter(rp.mesh1Domain + "_" + currentIsoDate.getYear() + String.format("%02d", currentIsoDate.getMonth()) + String.format("%02d", currentIsoDate.getDay()) + "*_rev.nc"),
//                                    null);
//                            ISO_datestr yesterday = ISO_datestr.getYesterday(currentIsoDate);
//                            if (isLastDay) {
//                                System.out.println("** Last day - reading same hydro file twice **");
//                                yesterday = currentIsoDate;
//                            }
//                            List<File> files2 = (List<File>) FileUtils.listFiles(
//                                    new File(rp.datadir + rp.datadirPrefix + yesterday.getYear() + rp.datadirSuffix + FileSystems.getDefault().getSeparator()),
//                                    new WildcardFileFilter(rp.mesh1Domain + "_" + yesterday.getYear() + String.format("%02d", yesterday.getMonth()) + String.format("%02d", yesterday.getDay()) + "*_rev.nc"),
//                                    null);
//                            // Read both files and combine
//                            hydroFields.add(new HydroField(files1.get(0).getCanonicalPath(), files2.get(0).getCanonicalPath(), varNames1, null, null, null, "FVCOM", rp));
//                        }
//                        readAttempt = 5; // success, so exit loop
//                    } catch (Exception ignored) {
//                        System.out.printf("Failed reading hydrofile on attempt %d. Retrying...\n", readAttempt+1);
//                        readAttempt++;
//                        Thread.sleep(5000);
//                        if (readAttempt == 5) {
//                            System.out.println("Hydro file not found, check PROPERTIES: datadir, datadirPrefix, datadirSuffix, location, minchVersion");
//                            if (!rp.backwards) {
//                                System.err.println("Requested file: " + rp.datadir + rp.datadirPrefix + currentIsoDate.getYear() + rp.datadirSuffix + FileSystems.getDefault().getSeparator()
//                                        + rp.mesh1Domain + "_" + currentIsoDate.getYear() + String.format("%02d", currentIsoDate.getMonth()) + String.format("%02d", currentIsoDate.getDay()) + "*.nc");
//                            } else {
//                                System.err.println("Requested file: " + rp.datadir + rp.datadirPrefix + currentIsoDate.getYear() + rp.datadirSuffix + FileSystems.getDefault().getSeparator()
//                                        + rp.mesh1Domain + "_" + currentIsoDate.getYear() + String.format("%02d", currentIsoDate.getMonth()) + String.format("%02d", currentIsoDate.getDay()) + "*_rev.nc");
//                            }
//                            System.exit(1);
//                        }
//                    }
//                }
//            } else if (mesh.getType().equalsIgnoreCase("ROMS_TRI")) {
//                String filename1 = rp.datadir2 + rp.datadir2Prefix + currentIsoDate.getYear() + rp.datadir2Suffix + FileSystems.getDefault().getSeparator()
//                        + "NEATL_" + currentIsoDate.getYear() + String.format("%02d", currentIsoDate.getMonth()) + String.format("%02d", currentIsoDate.getDay()) + ".nc";
//                String filename2 = rp.datadir2 + rp.datadir2Prefix + currentIsoDate.getYear() + rp.datadir2Suffix + FileSystems.getDefault().getSeparator()
//                        + "NEATL_" + currentIsoDate.getYear() + String.format("%02d", currentIsoDate.getMonth()) + String.format("%02d", currentIsoDate.getDay()) + ".nc";
//                String[] varNames1 = {"u", "v", "", "", ""};
//                // Read both files and combine
//                hydroFields.add(new HydroField(filename1, filename2, varNames1, null, null, null, "ROMS_TRI", rp));
//            }
//
//            m++;
//
//        }
//        return hydroFields;
//    }


    public static HydroField readHydroField(Mesh mesh, ISO_datestr currentIsoDate, RunProperties rp, int m) throws InterruptedException {
        if (mesh.getType().equalsIgnoreCase("FVCOM")) {
            // Occasionally the read fails for unknown reasons when the .nc file exists, so try up to 5 times before quitting
            int readAttempt = 0;
            while (readAttempt < 5) {
                try {
                    String[] varNames1 = new String[]{"u", "v", "ww", "salinity", "temp", "zeta", "kh", "viscofh", "short_wave", "Hsig", "Dir", "Tm01"};

                    // Inelegant, but it works and I'm in a hurry. Clean it up if you're procrastinating on something else.
                    String[] dirsT1 = new String[]{
                            rp.hfDir0 + rp.hfDirPrefix0 + currentIsoDate.getYear() + rp.hfDirSuffix0 + FileSystems.getDefault().getSeparator(),
                            rp.hfDir1 + rp.hfDirPrefix1 + currentIsoDate.getYear() + rp.hfDirSuffix1 + FileSystems.getDefault().getSeparator(),
                            rp.hfDir2 + rp.hfDirPrefix2 + currentIsoDate.getYear() + rp.hfDirSuffix2 + FileSystems.getDefault().getSeparator()};
                    String[] filesT1 = new String[]{
                            rp.hfFilePrefix0 + "_" + currentIsoDate.getYear() + String.format("%02d", currentIsoDate.getMonth()) + String.format("%02d", currentIsoDate.getDay()) + "*.nc",
                            rp.hfFilePrefix1 + "_" + currentIsoDate.getYear() + String.format("%02d", currentIsoDate.getMonth()) + String.format("%02d", currentIsoDate.getDay()) + "*.nc",
                            rp.hfFilePrefix2 + "_" + currentIsoDate.getYear() + String.format("%02d", currentIsoDate.getMonth()) + String.format("%02d", currentIsoDate.getDay()) + "*.nc"};
                    List<File> files1 = (List<File>) FileUtils.listFiles(new File(dirsT1[m]), new WildcardFileFilter(filesT1[m]), null);
                    // Read both files and combine
                    if (m == 2) {
                        return(new HydroField(files1.get(0).getCanonicalPath(), varNames1, null, null, null, true, rp));
                    } else {
                        return(new HydroField(files1.get(0).getCanonicalPath(), varNames1, null, null, null, rp));
                    }
                } catch (Exception ignored) {
                    System.out.printf("Failed reading hydrofile on attempt %d. Retrying...\n", readAttempt + 1);
                    readAttempt++;
                    Thread.sleep(5000);
                    if (readAttempt == 5) {
                        System.out.println("Hydro file not found, check PROPERTIES: datadir, datadirPrefix, datadirSuffix, location, minchVersion");
                        if (!rp.backwards) {
                            System.err.println("Requested file: " + rp.hfDir0 + rp.hfDirPrefix0 + currentIsoDate.getYear() + rp.hfDirSuffix0 + FileSystems.getDefault().getSeparator()
                                    + rp.hfFilePrefix0 + "_" + currentIsoDate.getYear() + String.format("%02d", currentIsoDate.getMonth()) + String.format("%02d", currentIsoDate.getDay()) + "*.nc");
                        } else {
                            System.err.println("Requested file: " + rp.hfDir0 + rp.hfDirPrefix0 + currentIsoDate.getYear() + rp.hfDirSuffix0 + FileSystems.getDefault().getSeparator()
                                    + rp.hfFilePrefix0 + "_" + currentIsoDate.getYear() + String.format("%02d", currentIsoDate.getMonth()) + String.format("%02d", currentIsoDate.getDay()) + "*_rev.nc");
                        }
                        System.exit(1);
                    }
                }
            }
        } else if (mesh.getType().equalsIgnoreCase("ROMS_TRI")) {
//            String filename1 = rp.hfDir1 + rp.hfDirPrefix1 + currentIsoDate.getYear() + rp.hfDirSuffix1 + FileSystems.getDefault().getSeparator()
//                    + "NEATL_" + currentIsoDate.getYear() + String.format("%02d", currentIsoDate.getMonth()) + String.format("%02d", currentIsoDate.getDay()) + ".nc";
//            String[] varNames1 = {"u", "v", "", "", ""};
//            // Read both files and combine
//            return(new HydroField(filename1, varNames1, null, null, null, "ROMS_TRI", rp));
        }
        return(null);
    }

    // Note: Stokes drift variables are stored in a separate structure as relevant .nc files are assumed to be distinct
    public static HydroField mergeHydroFields(HydroField hf1, HydroField hf2, RunProperties rp) {
        // Get hydrodynamics
        float[][][] u1 = hf1.getU(), u2 = hf2.getU();
        float[][][] u = new float[u1.length + 1][u1[0].length][u1[0][0].length];
        float[][][] v1 = hf1.getV(), v2 = hf2.getV();
        float[][][] v = new float[v1.length + 1][v1[0].length][v1[0][0].length];
        float[][][] w1 = hf1.getW(), w2 = hf2.getW();
        float[][][] w = new float[w1.length + 1][w1[0].length][w1[0][0].length];

        int nHour = u1.length;
        int nDep = u1[0].length;
        int nElem = u1[0][0].length;
        int nNode = 0;

        float[][][] s = null, s1 = null, s2 = null, t = null, t1 = null, t2 = null, k = null, k1 = null, k2 = null, vh = null, vh1 = null, vh2 = null;
        float[][] zeta = null, zeta1 = null, zeta2 = null, light = null, light1 = null, light2 = null;
        if (rp.needS) {
            s1 = hf1.getS();
            s2 = hf2.getS();
            s = new float[s1.length + 1][s1[0].length][s1[0][1].length];
            nNode = s1[0][0].length;
        }
        if (rp.needT) {
            t1 = hf1.getT();
            t2 = hf2.getT();
            t = new float[t1.length + 1][t1[0].length][t1[0][0].length];
            nNode = t1[0][0].length;
        }
        if (rp.needZeta) {
            zeta1 = hf1.getZeta();
            zeta2 = hf2.getZeta();
            zeta = new float[zeta1.length + 1][zeta1[0].length];
            nNode = zeta1[0].length;
        }
        if (rp.needLight) {
            light1 = hf1.getLight();
            light2 = hf2.getLight();
            light = new float[light1.length + 1][light1[0].length];
        }
        if (rp.needK) {
            k1 = hf1.getK();
            k2 = hf2.getK();
            k = new float[k1.length + 1][k1[0].length][k1[0][0].length];
        }
        if (rp.needVh) {
            vh1 = hf1.getVh();
            vh2 = hf2.getVh();
            vh = new float[vh1.length + 1][vh1[0].length][vh1[0][0].length];
        }

        for (int hour = 0; hour < nHour; hour++) {
            for (int dep = 0; dep < nDep; dep++) {
                for (int elem = 0; elem < nElem; elem++) {
                    u[hour][dep][elem] = u1[hour][dep][elem];
                    v[hour][dep][elem] = v1[hour][dep][elem];
                    w[hour][dep][elem] = w1[hour][dep][elem];
                }
                for (int node = 0; node < nNode; node++) {
                    if (s != null) {
                        s[hour][dep][node] = s1[hour][dep][node];
                    }
                    if (t != null) {
                        t[hour][dep][node] = t1[hour][dep][node];
                    }
                    if (dep == 0) {
                        if (zeta != null) {
                            zeta[hour][node] = zeta1[hour][node];
                        }
                        if (light != null) {
                            light[hour][node] = light1[hour][node];
                        }
                    }
                    if (k != null) {
                        k[hour][dep][node] = k1[hour][dep][node];
                        // easiest way to adjust for extra sigma level
                        if (dep == u1[0].length - 1) {
                            k[hour][dep + 1][node] = k1[hour][dep + 1][node];
                        }
                    }
                    if (vh != null) {
                        vh[hour][dep][node] = vh1[hour][dep][node];
                    }
                }
            }
        }
        for (int dep = 0; dep < nDep; dep++) {
            for (int elem = 0; elem < nElem; elem++) {
                u[u.length - 1][dep][elem] = u2[0][dep][elem];
                v[u.length - 1][dep][elem] = v2[0][dep][elem];
                w[u.length - 1][dep][elem] = w2[0][dep][elem];
            }
            for (int node = 0; node < nNode; node++) {
                if (s != null) {
                    s[u.length - 1][dep][node] = s2[0][dep][node];
                }
                if (t != null) {
                    t[u.length - 1][dep][node] = t2[0][dep][node];
                }
                if (dep == 0) {
                    if (zeta != null) {
                        zeta[u.length - 1][node] = zeta2[0][node];
                    }
                    if (light != null) {
                        light[u.length - 1][node] = light2[0][node];
                    }
                }
                if (k != null) {
                    k[u.length - 1][dep][node] = k2[0][dep][node];
                    // easiest way to adjust for extra sigma level
                    if (dep == u1[0].length - 1) {
                        k[u.length - 1][dep + 1][node] = k2[0][dep + 1][node];
                    }
                }
                if (vh != null) {
                    vh[u.length - 1][dep][node] = vh2[0][dep][node];
                }
            }
        }
        return(new HydroField(u, v, w, s, t, zeta, k, vh, light, null, null, null));
    }


    public static HydroField mergeHydroFieldsStokes(HydroField hf1, HydroField hf2, RunProperties rp) {
        // Get hydrodynamics
        float[][] Hsig1 = hf1.getHsig(), Hsig2 = hf2.getHsig();
        float[][] Hsig = new float[Hsig1.length + 1][Hsig1[0].length];
        float[][] Dir1 = hf1.getDir(), Dir2 = hf2.getDir();
        float[][] Dir = new float[Dir1.length + 1][Dir1[0].length];
        float[][] Tm1 = hf1.getTm(), Tm2 = hf2.getTm();
        float[][] Tm = new float[Tm1.length + 1][Tm1[0].length];

        int nHour = Hsig1.length;
        int nNode = Hsig1[0].length;

        for (int hour = 0; hour < nHour; hour++) {
            for (int node = 0; node < nNode; node++) {
                Hsig[hour][node] = Hsig1[hour][node];
                Dir[hour][node] = Dir1[hour][node];
                Tm[hour][node] = Tm1[hour][node];
            }
        }
        for (int node = 0; node < nNode; node++) {
            Hsig[Hsig.length - 1][node] = Hsig2[0][node];
            Dir[Dir.length - 1][node] = Dir2[0][node];
            Tm[Tm.length - 1][node] = Tm2[0][node];
        }
        return (new HydroField(null, null, null, null, null, null, null, null, null, Hsig, Dir, Tm));
    }


    /**
     * Count the number of particles in different states (free, viable, settled,
     * exited domain)
     */
    public static int[] particleCounts(List<Particle> parts) {
        int[] freeViableSettleExit = new int[4];
        for (Particle p : parts) {
            freeViableSettleExit[0] += p.isFree() ? 1 : 0;
            freeViableSettleExit[1] += p.isInfectious() ? 1 : 0;
            freeViableSettleExit[2] += p.hasArrived() ? 1 : 0;
            freeViableSettleExit[3] += p.hasExited() ? 1 : 0;
        }
        return freeViableSettleExit;
    }

    public static int[] particleCountStatus(List<Particle> parts, int[] previousStatus) {
        int[] status = new int[6];
        for (Particle p : parts) {
            int index = switch(p.getStatus()) {
                case 0 -> 0;
                case 1 -> 1;
                case 2 -> 2;
                case 3 -> 3;
                case 66 -> 4;
                case 666 -> 5;
                default -> 10;
            };
            status[index] += 1;
        }
        status[3] += previousStatus[3];
        status[4] += previousStatus[4];
        status[5] += previousStatus[5];
        return status;
    }


    public static int[] nonZeroVals(int[] A) {
        int count = 0;
        for (int j : A) {
            if (j != 0) {
                count++;
            }
        }
        int[] temp = new int[count];
        int p = 0;
        for (int j : A) {
            if (j != 0) {
                temp[p++] = j;
            }
        }
        return temp;
    }

    public static float[][] nonZeroRows(float[][] A) {
        int count = 0;
        List<Integer> list = new ArrayList<>();

        for (int i = 0; i < A.length; i++) {
            // Check whether ANY of the elements on this row !=0
            for (int j = 0; j < A[0].length; j++) {
                if (A[i][j] > 0) {
                    list.add(i);
                    count++;
                    break;
                }
            }
        }
        float[][] temp = null;
        if (count > 0) {
            temp = new float[count][A[0].length + 1];
            int p = 0;
            for (int row : list) {
                temp[p][0] = row;
                for (int j = 0; j < A[0].length; j++) {
                    if (A[row][j] > 0) {
                        temp[p][j + 1] = A[row][j];
                    }
                }
                p++;
            }
        }
        return temp;
    }


    public static void memTest() {
        long heapSize = Runtime.getRuntime().totalMemory();
        System.out.println("Total heap memory " + heapSize);
        long heapFreeSize = Runtime.getRuntime().freeMemory();
        System.out.println("Free heap memory " + heapFreeSize);

    }

    public void setupOutput() {
    }

    public void writeOutput() {
    }
}
