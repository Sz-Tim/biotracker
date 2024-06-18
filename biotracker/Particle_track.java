package biotracker;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.filefilter.WildcardFileFilter;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.nio.file.FileSystems;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Date;
import java.util.List;
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
        System.out.println("-----------------------------------------------------------");
        System.out.printf("Location           = %s\n", rp.location);
        System.out.printf("N_parts/site       = %d\n", rp.nparts);
        System.out.printf("hydromod dt (s)    = %.3f\n", rp.dt);
        System.out.printf("hydromod rec/file  = %d\n", rp.recordsPerFile1);
        System.out.printf("stepsperstep       = %d\n", rp.stepsPerStep);
        System.out.printf("firstfile          = %s\n", rp.start_ymd);
        System.out.printf("lastfile           = %s\n", rp.end_ymd);
        System.out.printf("Simulated dur. (d) = %d\n", numberOfDays);
        System.out.printf("Simulated dur. (s) = %d\n", numberOfDays * 86400);
        System.out.printf("RK4                = %s\n", rp.rk4);
        System.out.printf("Fixed depth        = %b\n", rp.fixDepth);
        System.out.printf("Max particle depth = %.3f\n", rp.maxDepth);
        System.out.printf("Viable time (h)    = %.3f\n", rp.viabletime);
        System.out.printf("Threshold distance = %d\n", rp.thresh);
        System.out.printf("Diffusion D_h      = %s\n", rp.variableDh ? "variable" : "" + rp.D_h);
        System.out.printf("Diffusion D_v      = %s\n", rp.variableDhV ? "variable" : "" + rp.D_hVert);
        System.out.printf("EPSG:27700         = %b\n", rp.coordOS);
        System.out.println("-----------------------------------------------------------");

        // --------------------------------------------------------------------------------------
        // File reading and domain configuration
        // --------------------------------------------------------------------------------------       
        List<Mesh> meshes = new ArrayList<>();
        meshes.add(new Mesh(rp.mesh1, rp.mesh1Type, rp.coordOS));
        if (!rp.mesh2.isEmpty()) {
            meshes.add(new Mesh(rp.mesh2, rp.mesh2Type, rp.coordOS));
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
        List<HydroField> hydroFields = new ArrayList<>();
        List<HydroField> hfToday;
        List<HydroField> hfTomorrow = null;
        ArrayList<String> missingHydroFiles = IOUtils.checkHydroFilesExist(rp, currentIsoDate, endIsoDate, numberOfDays);
        if(!missingHydroFiles.isEmpty()) {
            System.out.println("\nWarning! Cannot find the following hydrodynamic files:");
            for (String missing: missingHydroFiles) {
                System.out.println(missing);
            }
            System.out.println("-- The previous day will be used for each instead!");
        } else {
            System.out.println("  All found");
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
        // Final setup bits
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

        int stepcount = 0;
        double elapsedHours = 0; // updated in HOURS as the simulation progresses

        int[] freeViableSettleExit;

        int numberOfExecutorThreads = rp.parallelThreads;
        System.out.println("Number of executor threads = " + numberOfExecutorThreads);
        ExecutorService executorService = Executors.newFixedThreadPool(numberOfExecutorThreads);
        CompletionService<List<Particle>> executorCompletionService = new ExecutorCompletionService<>(executorService);

        final Collection<Callable<List<Particle>>> callables = new ArrayList<>();

        String particleRestartHeader = "hour,ID,startDate,age,startLocation,x,y,elem,status,density,mesh,depth,depthLayer,degreeDays,xTot,yTot,xyTot,zTot";
        String arrivalHeader = "ID,startDate,startTime,startLocation,endDate,endTime,endLocation,age,density";

        // Set up arrays to hold particle density*hour counts
        int pstepsInd2 = 1;
        if (rp.splitPsteps) {
            pstepsInd2 = habitat.size();
        }
        float[][] pstepsImmature = new float[meshes.get(0).getNElems()][pstepsInd2];
        float[][] pstepsMature = new float[meshes.get(0).getNElems()][pstepsInd2];

        // Set up array to hold connectivity counts
        float[][] connectivity = new float[habitat.size()][habitatEnd.size()];

        int[][] elemActivity = new int[meshes.get(0).getNElems()][3]; // count of sink, swim, float within each element
        int[][] hourActivity = new int[numberOfDays*24+1][3]; // count of sink, swim, float within each hour
        float[][] vertDistrImmature = new float[meshes.get(0).getNElems()][rp.vertDistrMax+1]; // bin upper z limits: from 0 to vertDistrMax
        float[][] vertDistrMature = new float[meshes.get(0).getNElems()][rp.vertDistrMax+1]; // bin upper z limits: from 0 to vertDistrMax
        if (rp.recordMovement) {
            IOUtils.writeMovementsHeader("ID,date,hour,step,startDate,age,density,x,y,z,layer,status,degreeDays,sink,swim,temp,salinity,mortality,tempSurface,dX,dY,dZ",
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
                System.out.printf("\n------ Day %d (%s) - Stepcount %d (%.1f hrs) ------ \n",
                        fnum + 1, currentIsoDate, stepcount, elapsedHours);
                System.out.println("Elapsed time (s) = " + (splitTime - startTime) / 1000.0);
                System.out.println("Last 24hr time (s) = " + (splitTime - currTime) / 1000.0);
                currTime = System.currentTimeMillis();

                // COUNT the number of particles in different states
                freeViableSettleExit = particleCounts(particles);
                System.out.println("Free particles    = " + freeViableSettleExit[0]);
                System.out.println("Viable particles  = " + freeViableSettleExit[1]);
                System.out.println("Arrival count     = " + freeViableSettleExit[2]);
                System.out.println("Boundary exits    = " + freeViableSettleExit[3]);


                  // default, run loop forwards
                // ---- LOOP OVER ENTRIES IN THE HYDRO OUTPUT ------------------------
                for (int currentHour = 0; currentHour < 24; currentHour++) {

                    System.out.printf("--------- Day %d, %dh ----------\n", fnum+1, currentHour);
                    // Calculate current time of the day (complete hours elapsed since midnight)
                    if (!rp.daylightPath.isEmpty()) {
                        isDaytime = daylightHours[fnum][0] <= currentHour && daylightHours[fnum][1] >= currentHour;
                    }

                    // Read new hydrodynamic fields?
                    if (currentHour == 0) {
                        String currentDateShort = currentIsoDate.toStringShort();
                        ISO_datestr nextIsoDate = ISO_datestr.getTomorrow(currentIsoDate);
                        String nextDateShort = nextIsoDate.toStringShort();
                        boolean missingToday = missingHydroFiles.contains(currentDateShort);
                        boolean missingTomorrow = missingHydroFiles.contains(nextDateShort);
                        if (missingToday) {
                            System.out.print("Can't find file: Reusing hydrodynamics from yesterday\n");
                        } else {
                            hfToday = fnum == 0 ? readHydroField(meshes, currentIsoDate, rp) : hfTomorrow;
                            hfTomorrow = (missingTomorrow | isLastDay) ? hfToday : readHydroField(meshes, nextIsoDate, rp);
                            hydroFields = mergeHydroFields(hfToday, hfTomorrow, rp);
                        }
                    }

                    for (int i=0; i < habitat.size(); i++) {
                        int siteElem = habitat.get(i).getContainingFVCOMElem();
                        int m = habitat.get(i).getContainingMesh();
                        int nLayers = (int) Mesh.findNearestSigmas(30.0, meshes.get(m).getSiglay(), (float) habitat.get(i).getDepth())[0][0];
                        double[][] currentConditions = new double[nLayers][6];
                        double[] siteLoc = new double[2];
                        siteLoc[0] = habitat.get(i).getLocation()[0];
                        siteLoc[1] = habitat.get(i).getLocation()[1];
                        for (int j = 0; j < nLayers; j++) {
                            currentConditions[j][0] = hydroFields.get(m).getU()[currentHour][j][siteElem];
                            currentConditions[j][1] = hydroFields.get(m).getV()[currentHour][j][siteElem];
                            currentConditions[j][2] = hydroFields.get(m).getW()[currentHour][j][siteElem];
                            currentConditions[j][3] = Math.sqrt(currentConditions[j][0]*currentConditions[j][0] + currentConditions[j][1]*currentConditions[j][1]);
                            currentConditions[j][4] = rp.needS ? hydroFields.get(m).getAvgFromTrinodes(meshes.get(m), siteLoc, j, siteElem, currentHour, "salinity", rp) : -9999;
                            currentConditions[j][5] = rp.needT ? hydroFields.get(m).getAvgFromTrinodes(meshes.get(m), siteLoc, j, siteElem, currentHour, "temp", rp) : -9999;
                            habitat.get(i).addEnvCondition(currentConditions[j]);
                        }
                        // calculate eggs per female per timestep based on temperature and scale initial particle densities
                        // Note: Norwegian model uses Stien et al 2005: N_naup = N_fish * N_female * 0.17 * (temp+4.28)^2
                        // Linear model comes from Kragesteen
                        // These are nearly identical under the temperatures in Scotland (r = 0.998)
                        float[][] eggSigmas = Mesh.findNearestSigmas(2, meshes.get(m).getSiglay(), (float) habitat.get(i).getDepth());
                        double eggTemperature = hydroFields.get(m).getValueAtDepth(meshes.get(m), siteElem, siteLoc, 2, currentHour, "temp", rp, eggSigmas);
                        double eggsPerFemale = rp.eggTemp_b0 * (rp.dt / (60*60*24));
                        if (rp.eggTemp_b1 > 0.001) {
                            // eggsPerFemale = (rp.eggTemp_b0 + rp.eggTemp_b1 * eggTemperature) / (60*60*24/rp.dt);
                            eggsPerFemale = rp.eggTemp_b0 * Math.pow((eggTemperature + rp.eggTemp_b1), 2) * (rp.dt / (60*60*24));
                        }
                        habitat.get(i).setScale(siteDensities[i][fnum] / rp.nparts * eggsPerFemale);
                    }

                    // Create new particles, if releases are scheduled hourly, or if release is scheduled for this hour
                    double hourRemainder = Math.round((elapsedHours % rp.releaseInterval) * 100.0) / 100.0;
                    if ((rp.releaseScenario == 0 && elapsedHours >= rp.releaseTime && allowRelease) ||
                            rp.releaseScenario == 1 && (hourRemainder == 0 || hourRemainder == rp.releaseInterval) ||
                            (rp.releaseScenario == 2 && elapsedHours >= rp.releaseTime && elapsedHours <= rp.releaseTimeEnd)) {
                        List<Particle> newParts = createNewParticles(habitat, meshes, rp, currentIsoDate, currentHour, numParticlesCreated);
                        particles.addAll(newParts);
                        numParticlesCreated = numParticlesCreated + (rp.nparts * habitat.size());
                        System.out.printf("  %,d new particles (%,d active of %,d total)\n", newParts.size(), particles.size(), numParticlesCreated);
                        // If only one release to be made, prevent further releases
                        if (rp.releaseScenario == 0) {
                            allowRelease = false;
                        }
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

                    if (rp.recordLocations) {
                        IOUtils.particlesToRestartFile(particles, currentHour, "locations_" + today + ".csv", true, rp, 1); // rp.nparts * rp.numberOfDays * 10
                    }


                    double elapsedHourRemainder = Math.round((elapsedHours % 1) * 100.0) / 100.0;
                    if (elapsedHourRemainder == 1 || elapsedHourRemainder == 0) {
                        // It's the end of an hour, so if particles are allowed to infect more than once, reactivate them
                        for (Particle part : particles) {
                            if (part.hasSettledThisHour()) {
                                // Save arrival
                                if (rp.recordArrivals) {
                                    IOUtils.arrivalToFile(part, currentIsoDate, currentHour, "arrivals_" + today + ".csv", true);
                                }
                                // Add arrival to connectivity file
                                int destIndex = siteEndNames.indexOf(part.getLastArrival());
                                connectivity[part.getStartIndex()][destIndex] += (float) part.getDensity();

                                // Reset ability to settle
                                part.setSettledThisHour(false);
                            }
                        }
                    }


                    // Hourly updates to pstep arrays
                    if (rp.recordPsteps) {
                        IOUtils.pstepsUpdater(particles, rp, pstepsMature, pstepsImmature, rp.dt);
                    }

                    // Write Psteps
                    if (rp.recordPsteps && stepcount % (rp.pstepsInterval * rp.stepsPerStep) == 0) {
                        System.out.println("Writing psteps");
                        if (rp.recordImmature) {
                            IOUtils.writeNonZeros2DArrayToCSV(pstepsImmature, "i,NA,value","%d,%d,%.4e" , "pstepsImmature_" + today + "_" + Math.round(elapsedHours) + ".csv");
                            pstepsImmature = new float[meshes.get(0).getNElems()][pstepsInd2];
                        }
                        IOUtils.writeNonZeros2DArrayToCSV(pstepsMature, "i,NA,value","%d,%d,%.4e", "pstepsMature_" + today + "_" + Math.round(elapsedHours) + ".csv");
                        pstepsMature = new float[meshes.get(0).getNElems()][pstepsInd2];
                    }

                    if (rp.recordConnectivity && stepcount % (rp.connectivityInterval * rp.stepsPerStep) == 0) {
                        IOUtils.writeNonZeros2DArrayToCSV(connectivity, "source,destination,value", "%d,%d,%.4e", "connectivity_" + today + "_" + Math.round(elapsedHours) + ".csv");
                        connectivity = new float[habitat.size()][habitatEnd.size()];
                    }

                    if (rp.recordVertDistr) {
                        IOUtils.vertDistrUpdater(particles, rp, vertDistrMature, vertDistrImmature, rp.dt);
                        if (stepcount % (rp.vertDistrInterval * rp.stepsPerStep) == 0) {
                            System.out.println("Writing vertical distribution");
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
        out.println("site,x,y,depth,mesh,centroid,elem,MeshType,u,v,w,uv,salinity,temperature,k");
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
                    rp.mortalityRate, startDensity, currentDate, currentTime, rp.coordOS, rp.species);
            p.setMesh(meshStart);
            p.setElem(elemFVCOMStart);
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
    public static List<HydroField> readHydroFields(List<Mesh> meshes, ISO_datestr currentIsoDate, int currentHour, boolean isLastDay, RunProperties rp) throws InterruptedException {
        List<HydroField> hydroFields = new ArrayList<>();

        // 24 hr files only case - read once a day
        int m = 0;
        for (Mesh mesh : meshes) {
            if (mesh.getType().equalsIgnoreCase("FVCOM")) {

                // Occasionally the read fails for unknown reasons when the .nc file exists, so try up to 5 times before quitting
                int readAttempt = 0;
                while (readAttempt < 5) {
                    try{
                        // Dima file naming format: minch2_20171229_0003.nc
                        String[] varNames1 = new String[]{"u", "v", "ww", "salinity", "temp", "zeta", "kh", "viscofh", "short_wave"};

                        // Normal "forwards time"
                        if (!rp.backwards) {
                            ISO_datestr tomorrow = ISO_datestr.getTomorrow(currentIsoDate);
                            if (isLastDay) {
                                System.out.println("** Last day or missing day - reading same hydro file twice **");
                                tomorrow = currentIsoDate;
                            }
                            // Inelegant, but it works and I'm in a hurry. Clean it up if you're procrastinating on something else.
                            String[] dirsT1 = new String[]{
                                    rp.datadir + rp.datadirPrefix + currentIsoDate.getYear() + rp.datadirSuffix + FileSystems.getDefault().getSeparator(),
                                    rp.datadir2 + rp.datadir2Prefix + currentIsoDate.getYear() + rp.datadir2Suffix + FileSystems.getDefault().getSeparator()};
                            String[] dirsT2 = new String[]{
                                    rp.datadir + rp.datadirPrefix + tomorrow.getYear() + rp.datadirSuffix + FileSystems.getDefault().getSeparator(),
                                    rp.datadir2 + rp.datadir2Prefix + tomorrow.getYear() + rp.datadir2Suffix + FileSystems.getDefault().getSeparator()};
                            String[] filesT1 = new String[]{
                                    rp.location + rp.minchVersion + "_" + currentIsoDate.getYear() + String.format("%02d", currentIsoDate.getMonth()) + String.format("%02d", currentIsoDate.getDay()) + "*.nc",
                                    rp.location2 + rp.minchVersion2 + "_" + currentIsoDate.getYear() + String.format("%02d", currentIsoDate.getMonth()) + String.format("%02d", currentIsoDate.getDay()) + "*.nc"};
                            String[] filesT2 = new String[]{
                                    rp.location + rp.minchVersion + "_" + tomorrow.getYear() + String.format("%02d", tomorrow.getMonth()) + String.format("%02d", tomorrow.getDay()) + "*.nc",
                                    rp.location2 + rp.minchVersion2 + "_" + tomorrow.getYear() + String.format("%02d", tomorrow.getMonth()) + String.format("%02d", tomorrow.getDay()) + "*.nc"};
                            List<File> files1 = (List<File>) FileUtils.listFiles(new File(dirsT1[m]), new WildcardFileFilter(filesT1[m]), null);
                            List<File> files2 = (List<File>) FileUtils.listFiles(new File(dirsT2[m]), new WildcardFileFilter(filesT2[m]), null);
                            // Read both files and combine
                            hydroFields.add(new HydroField(files1.get(0).getCanonicalPath(), files2.get(0).getCanonicalPath(), varNames1, null, null, null, "FVCOM", rp));
                        }
                        // Instead read time backwards, so need yesterday instead
                        else {
                            List<File> files1 = (List<File>) FileUtils.listFiles(
                                    new File(rp.datadir + rp.datadirPrefix + currentIsoDate.getYear() + rp.datadirSuffix + FileSystems.getDefault().getSeparator()),
                                    new WildcardFileFilter(rp.location + rp.minchVersion + "_" + currentIsoDate.getYear() + String.format("%02d", currentIsoDate.getMonth()) + String.format("%02d", currentIsoDate.getDay()) + "*_rev.nc"),
                                    null);
                            ISO_datestr yesterday = ISO_datestr.getYesterday(currentIsoDate);
                            if (isLastDay) {
                                System.out.println("** Last day - reading same hydro file twice **");
                                yesterday = currentIsoDate;
                            }
                            List<File> files2 = (List<File>) FileUtils.listFiles(
                                    new File(rp.datadir + rp.datadirPrefix + yesterday.getYear() + rp.datadirSuffix + FileSystems.getDefault().getSeparator()),
                                    new WildcardFileFilter(rp.location + rp.minchVersion + "_" + yesterday.getYear() + String.format("%02d", yesterday.getMonth()) + String.format("%02d", yesterday.getDay()) + "*_rev.nc"),
                                    null);
                            // Read both files and combine
                            hydroFields.add(new HydroField(files1.get(0).getCanonicalPath(), files2.get(0).getCanonicalPath(), varNames1, null, null, null, "FVCOM", rp));
                        }
                        readAttempt = 5; // success, so exit loop
                    } catch (Exception ignored) {
                        System.out.printf("Failed reading hydrofile on attempt %d. Retrying...\n", readAttempt+1);
                        readAttempt++;
                        Thread.sleep(5000);
                        if (readAttempt == 5) {
                            System.out.println("Hydro file not found, check PROPERTIES: datadir, datadirPrefix, datadirSuffix, location, minchVersion");
                            if (!rp.backwards) {
                                System.err.println("Requested file: " + rp.datadir + rp.datadirPrefix + currentIsoDate.getYear() + rp.datadirSuffix + FileSystems.getDefault().getSeparator()
                                        + rp.location + rp.minchVersion + "_" + currentIsoDate.getYear() + String.format("%02d", currentIsoDate.getMonth()) + String.format("%02d", currentIsoDate.getDay()) + "*.nc");
                            } else {
                                System.err.println("Requested file: " + rp.datadir + rp.datadirPrefix + currentIsoDate.getYear() + rp.datadirSuffix + FileSystems.getDefault().getSeparator()
                                        + rp.location + rp.minchVersion + "_" + currentIsoDate.getYear() + String.format("%02d", currentIsoDate.getMonth()) + String.format("%02d", currentIsoDate.getDay()) + "*_rev.nc");
                            }
                            System.exit(1);
                        }
                    }
                }
            } else if (mesh.getType().equalsIgnoreCase("ROMS_TRI")) {
                String filename1 = rp.datadir2 + rp.datadir2Prefix + currentIsoDate.getYear() + rp.datadir2Suffix + FileSystems.getDefault().getSeparator()
                        + "NEATL_" + currentIsoDate.getYear() + String.format("%02d", currentIsoDate.getMonth()) + String.format("%02d", currentIsoDate.getDay()) + ".nc";
                String filename2 = rp.datadir2 + rp.datadir2Prefix + currentIsoDate.getYear() + rp.datadir2Suffix + FileSystems.getDefault().getSeparator()
                        + "NEATL_" + currentIsoDate.getYear() + String.format("%02d", currentIsoDate.getMonth()) + String.format("%02d", currentIsoDate.getDay()) + ".nc";
                String[] varNames1 = {"u", "v", "", "", ""};
                // Read both files and combine
                hydroFields.add(new HydroField(filename1, filename2, varNames1, null, null, null, "ROMS_TRI", rp));
            }

            m++;

        }
        return hydroFields;
    }


    public static List<HydroField> readHydroField(List<Mesh> meshes, ISO_datestr currentIsoDate, RunProperties rp) throws InterruptedException {
        List<HydroField> hydroFields = new ArrayList<>();
        // 24 hr files only case - read once a day
        int m = 0;
        for (Mesh mesh : meshes) {
            if (mesh.getType().equalsIgnoreCase("FVCOM")) {

                // Occasionally the read fails for unknown reasons when the .nc file exists, so try up to 5 times before quitting
                int readAttempt = 0;
                while (readAttempt < 5) {
                    try{
                        // Dima file naming format: minch2_20171229_0003.nc
                        String[] varNames1 = new String[]{"u", "v", "ww", "salinity", "temp", "zeta", "kh", "viscofh", "short_wave"};

                        // Inelegant, but it works and I'm in a hurry. Clean it up if you're procrastinating on something else.
                        String[] dirsT1 = new String[]{
                                rp.datadir + rp.datadirPrefix + currentIsoDate.getYear() + rp.datadirSuffix + FileSystems.getDefault().getSeparator(),
                                rp.datadir2 + rp.datadir2Prefix + currentIsoDate.getYear() + rp.datadir2Suffix + FileSystems.getDefault().getSeparator()};
                        String[] filesT1 = new String[]{
                                rp.location + rp.minchVersion + "_" + currentIsoDate.getYear() + String.format("%02d", currentIsoDate.getMonth()) + String.format("%02d", currentIsoDate.getDay()) + "*.nc",
                                rp.location2 + rp.minchVersion2 + "_" + currentIsoDate.getYear() + String.format("%02d", currentIsoDate.getMonth()) + String.format("%02d", currentIsoDate.getDay()) + "*.nc"};
                        List<File> files1 = (List<File>) FileUtils.listFiles(new File(dirsT1[m]), new WildcardFileFilter(filesT1[m]), null);
                        // Read both files and combine
                        hydroFields.add(new HydroField(files1.get(0).getCanonicalPath(), varNames1, null, null, null, "FVCOM", rp));
                        readAttempt = 5; // success, so exit loop
                    } catch (Exception ignored) {
                        System.out.printf("Failed reading hydrofile on attempt %d. Retrying...\n", readAttempt+1);
                        readAttempt++;
                        Thread.sleep(5000);
                        if (readAttempt == 5) {
                            System.out.println("Hydro file not found, check PROPERTIES: datadir, datadirPrefix, datadirSuffix, location, minchVersion");
                            if (!rp.backwards) {
                                System.err.println("Requested file: " + rp.datadir + rp.datadirPrefix + currentIsoDate.getYear() + rp.datadirSuffix + FileSystems.getDefault().getSeparator()
                                        + rp.location + rp.minchVersion + "_" + currentIsoDate.getYear() + String.format("%02d", currentIsoDate.getMonth()) + String.format("%02d", currentIsoDate.getDay()) + "*.nc");
                            } else {
                                System.err.println("Requested file: " + rp.datadir + rp.datadirPrefix + currentIsoDate.getYear() + rp.datadirSuffix + FileSystems.getDefault().getSeparator()
                                        + rp.location + rp.minchVersion + "_" + currentIsoDate.getYear() + String.format("%02d", currentIsoDate.getMonth()) + String.format("%02d", currentIsoDate.getDay()) + "*_rev.nc");
                            }
                            System.exit(1);
                        }
                    }
                }
            } else if (mesh.getType().equalsIgnoreCase("ROMS_TRI")) {
                String filename1 = rp.datadir2 + rp.datadir2Prefix + currentIsoDate.getYear() + rp.datadir2Suffix + FileSystems.getDefault().getSeparator()
                        + "NEATL_" + currentIsoDate.getYear() + String.format("%02d", currentIsoDate.getMonth()) + String.format("%02d", currentIsoDate.getDay()) + ".nc";
                String[] varNames1 = {"u", "v", "", "", ""};
                // Read both files and combine
                hydroFields.add(new HydroField(filename1, varNames1, null, null, null, "ROMS_TRI", rp));
            }
            m++;
        }
        return hydroFields;
    }

    public static List<HydroField> mergeHydroFields(List<HydroField> hf1, List<HydroField> hf2, RunProperties rp) {
        List<HydroField> hydroFields = new ArrayList<>();

        for (int m=0; m < hf1.size(); m++) {
            // Get hydrodynamics
            float[][][] u1 = hf1.get(m).getU(), u2 = hf2.get(m).getU();
            float[][][] u = new float[u1.length + 1][u1[0].length][u1[0][0].length];
            float[][][] v1 = hf1.get(m).getV(), v2 = hf2.get(m).getV();
            float[][][] v = new float[u1.length + 1][u1[0].length][u1[0][0].length];

            int nHour = u1.length;
            int nDep = u1[0].length;
            int nElem = u1[0][0].length;
            int nNode = 0;

            float[][][] w = null, w1 = null, w2 = null, s = null, s1 = null, s2 = null, t = null, t1 = null, t2 = null, k = null, k1 = null, k2 = null, vh = null, vh1 = null, vh2 = null;
            float[][] zeta = null, zeta1 = null, zeta2 = null, light = null, light1 = null, light2 = null;
            w1 = hf1.get(m).getW();
            w2 = hf2.get(m).getW();
            w = new float[w1.length + 1][w1[0].length][w1[0][0].length];
            if(rp.needS) {
                s1 = hf1.get(m).getS();
                s2 = hf2.get(m).getS();
                s = new float[s1.length + 1][s1[0].length][s1[0][1].length];
                nNode = s1[0][0].length;
            }
            if(rp.needT) {
                t1 = hf1.get(m).getT();
                t2 = hf2.get(m).getT();
                t = new float[t1.length + 1][t1[0].length][t1[0][0].length];
                nNode = t1[0][0].length;
            }
            if(rp.needZeta) {
                zeta1 = hf1.get(m).getZeta();
                zeta2 = hf2.get(m).getZeta();
                zeta = new float[zeta1.length + 1][zeta1[0].length];
                nNode = zeta1[0].length;
            }
            if(rp.needLight) {
                light1 = hf1.get(m).getLight();
                light2 = hf2.get(m).getLight();
                light = new float[light1.length + 1][light1[0].length];
            }
            if(rp.needK) {
                k1 = hf1.get(m).getK();
                k2 = hf2.get(m).getK();
                k = new float[k1.length + 1][k1[0].length][k1[0][0].length];
            }
            if(rp.needVh) {
                vh1 = hf1.get(m).getVh();
                vh2 = hf2.get(m).getVh();
                vh = new float[vh1.length + 1][vh1[0].length][vh1[0][0].length];
            }

            double sumU = 0;

            // Use this to test zero velocity/diffusion cases
            boolean createTest = false;
            float testU = 0;
            float testV = (float) 0.1;
            float testW = (float) 0.1;
            //noinspection ConstantConditions
            if (createTest) {
                for (int hour = 0; hour < u.length; hour++) {
                    for (int dep = 0; dep < u[0].length; dep++) {
                        for (int elem = 0; elem < u[0][0].length; elem++) {
                            u[hour][dep][elem] = testU;
                            v[hour][dep][elem] = testV;
                            w[hour][dep][elem] = testW;
                            sumU += u[hour][dep][elem];
                        }
                        if (!rp.readHydroVelocityOnly) {
                            for (int node = 0; node < zeta1[0].length; node++) {

                                if (s2 != null) {
                                    s[u.length - 1][dep][node] = 0;
                                }
                                if (t2 != null) {
                                    t[u.length - 1][dep][node] = 0;
                                }
                                if (k2 != null) {
                                    k[u.length - 1][dep][node] = 0;
                                    // easiest way to adjust for extra sigma level
                                    if (dep == u[0].length - 1) {
                                        k[u.length - 1][dep + 1][node] = 0;
                                    }
                                }
                                if (vh2 != null) {
                                    vh[u.length - 1][dep][node] = 0;
                                }
                                if (dep == 0) {
                                    zeta[u.length - 1][node] = 0;
                                    light[u.length - 1][node] = 0;
                                }
                            }
                        }
                    }
                }
            } else {
                for (int hour = 0; hour < nHour; hour++) {
                    for (int dep = 0; dep < nDep; dep++) {
                        for (int elem = 0; elem < nElem; elem++) {
                            u[hour][dep][elem] = u1[hour][dep][elem];
                            v[hour][dep][elem] = v1[hour][dep][elem];
                            if (w != null) {
                                w[hour][dep][elem] = w1[hour][dep][elem];
                            }
                            sumU += u[hour][dep][elem];
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
                        if (w != null) {
                            w[u.length - 1][dep][elem] = w2[0][dep][elem];
                        }
                        sumU += u[u.length - 1][dep][elem];
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
            }
            hydroFields.add(new HydroField(u, v, w, s, t, zeta, k, vh, light));
        }
        return hydroFields;
    }

    /**
     * Count the number of particles in different states (free, viable, settled,
     * exited domain)
     */
    public static int[] particleCounts(List<Particle> parts) {
        int[] freeViableSettleExit = new int[4];
        for (Particle p : parts) {
            freeViableSettleExit[0] += p.isFree() ? 1 : 0;
            freeViableSettleExit[1] += p.isViable() ? 1 : 0;
            freeViableSettleExit[2] += p.hasArrived() ? 1 : 0;
            freeViableSettleExit[3] += p.hasExited() ? 1 : 0;
        }
        return freeViableSettleExit;
    }


    // calculate a connectivity matrix detailing the 
    public static double[][] connectFromParticleArrivals(List<Particle> particles, int nStartLocs, int nEndLocs, int npartsPerSite) {
        double[][] connectMatrix = new double[nStartLocs][nEndLocs];
        for (Particle part : particles) {
            for (Arrival arrival : part.getArrivals()) {
                connectMatrix[arrival.getSourceLocation()][arrival.getArrivalLocation()] += arrival.getArrivalDensity() / npartsPerSite;
            }
        }
        return connectMatrix;
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
