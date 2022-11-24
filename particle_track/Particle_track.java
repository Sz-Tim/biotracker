package particle_track;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.filefilter.WildcardFileFilter;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
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

        System.out.println("Starting particle tracking program\n");
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
        System.out.printf("Starting densities = %s\n", rp.seasonalDensityPath.isEmpty() ? "constant" : "seasonal");
        System.out.printf("Vertical dynamics  = %b\n", rp.verticalDynamics);
        System.out.printf("Max particle depth = %.3f\n", rp.maxDepth);
        System.out.printf("Viable time (h)    = %.3f\n", rp.viabletime);
        System.out.printf("Threshold distance = %d\n", rp.thresh);
        System.out.printf("Diffusion D_h      = %.3f (diffusion: %s)\n", rp.D_h, rp.diffusion);
        System.out.printf("Diffusion D_v      = %s\n", rp.variableDiffusion ? "variable" : "" + rp.D_hVert);
        System.out.printf("Coord ref          = %s\n", rp.coordRef);
        System.out.println("-----------------------------------------------------------");

        // --------------------------------------------------------------------------------------
        // File reading and domain configuration
        // --------------------------------------------------------------------------------------       
        List<Mesh> meshes = new ArrayList<>();
        meshes.add(new Mesh(rp.mesh1, rp.mesh1Type, rp.coordRef));
        if (!rp.mesh2.equals("")) {
            meshes.add(new Mesh(rp.mesh2, rp.mesh2Type, rp.coordRef));
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
        ArrayList<String> missingHydroFiles = IOUtils.checkHydroFilesExist(rp, currentIsoDate, endIsoDate, numberOfDays);
        if(!missingHydroFiles.isEmpty()) {
            System.err.println("\nError! Cannot find the following hydrodynamic files:");
            for (String missing: missingHydroFiles) {
                System.err.println(missing);
            }
            System.exit(1);
        } else {
            System.out.println(": All found");
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
        double[] dailyDensities = new double[0];
        if (!rp.seasonalDensityPath.isEmpty()) {
            dailyDensities = IOUtils.readFileDouble1D(rp.seasonalDensityPath);
        }

        int stepcount = 0;
        double elapsedHours = 0; // updated in HOURS as the simulation progresses

        int[] freeViableSettleExit;

        int numberOfExecutorThreads = rp.parallelThreads;
        System.out.println("Number of executor threads = " + numberOfExecutorThreads);
        ExecutorService executorService = Executors.newFixedThreadPool(numberOfExecutorThreads);
        CompletionService<List<Particle>> executorCompletionService = new ExecutorCompletionService<>(executorService);

        final Collection<Callable<List<Particle>>> callables = new ArrayList<>();

        String particleRestartHeader = "hour ID startDate age startLocation x y elem status density mesh depth depthLayer degreeDays xTot yTot xyTot zTot";
        String arrivalHeader = "ID startDate startTime startLocation endDate endTime endLocation age density";

        // Set up arrays to hold particle density*hour counts
        int pstepsInd2 = 2;
        if (rp.splitPsteps) {
            pstepsInd2 = habitat.size();
        }
        float[][] pstepsImmature = new float[meshes.get(0).getNElems()][pstepsInd2];
        float[][] pstepsMature = new float[meshes.get(0).getNElems()][pstepsInd2];

        // Set up array to hold connectivity counts
        float[][] connectivity = new float[habitat.size()][habitatEnd.size()];

        int[][] elemActivity = new int[meshes.get(0).getNElems()][3]; // count of sink, swim, float within each element
        if (rp.recordMovement) {
            IOUtils.writeMovementsHeader("ID date hour step startDate age density x y z layer status degreeDays sink swim temp salinity mortality tempSurface dX dY dZ",
                    "movementFile.dat");
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
                    IOUtils.printFileHeader(particleRestartHeader, "locations_" + today + ".dat");
                }
                if (rp.recordArrivals) {
                    IOUtils.printFileHeader(arrivalHeader, "arrivals_" + today + ".dat");
                }

                if (!rp.seasonalDensityPath.isEmpty()) {
                    startDensity = dailyDensities[fnum];
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
                        hydroFields.clear();
                        hydroFields = readHydroFields(meshes, currentIsoDate, currentHour, isLastDay, rp);
                    }

                    for (HabitatSite site: habitat) {
                        int siteElem = site.getContainingFVCOMElem();
                        int nLayers = (int) Mesh.findNearestSigmas(30.0, meshes.get(0).getSiglay(), (float) site.getDepth())[0][0];
                        double[][] currentConditions = new double[nLayers][4];
                        double[] siteLoc = new double[2];
                        siteLoc[0] = site.getLocation()[0];
                        siteLoc[1] = site.getLocation()[1];
                        for (int i = 0; i < nLayers; i++) {
                            currentConditions[i][0] = hydroFields.get(0).getU()[currentHour][i][siteElem];
                            currentConditions[i][1] = hydroFields.get(0).getV()[currentHour][i][siteElem];
                            currentConditions[i][2] = hydroFields.get(0).getW()[currentHour][i][siteElem];
                            currentConditions[i][3] = hydroFields.get(0).getAvgFromTrinodes(meshes.get(0), siteLoc, i, siteElem, currentHour, "salinity", rp);
                            site.addEnvCondition(currentConditions[i]);
                        }
                    }

                    // Create new particles, if releases are scheduled hourly, or if release is scheduled for this hour
                    double hourRemainder = Math.round((elapsedHours % rp.releaseInterval) * 100.0) / 100.0;
                    if ((rp.releaseScenario == 0 && elapsedHours >= rp.releaseTime && allowRelease) ||
                            rp.releaseScenario == 1 && (hourRemainder == 0 || hourRemainder == rp.releaseInterval) ||
                            (rp.releaseScenario == 2 && elapsedHours >= rp.releaseTime && elapsedHours <= rp.releaseTimeEnd)) {
                        System.out.printf("  %dD movement, release density: %.3f\n", rp.verticalDynamics ? 3 : 2, startDensity);
                        List<Particle> newParts = createNewParticles(habitat, meshes, rp, currentIsoDate, currentHour, startDensity, numParticlesCreated);
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
                                        meshes, hydroFields, habitatEnd, allelems, isDaytime, today, elemActivity));
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
                                        meshes, hydroFields, habitatEnd, allelems, isDaytime, today, elemActivity);
                            }
                        }

                        elapsedHours += subStepDt / 3600.0;
                        stepcount++;
                    }

                    if (rp.recordLocations) {
                        IOUtils.particlesToRestartFile(particles, currentHour, "locations_" + today + ".dat", true, rp, 1); // rp.nparts * rp.numberOfDays * 10
                    }

                    // It's the end of an hour, so if particles are allowed to infect more than once, reactivate them
                    for (Particle part : particles) {
                        if (part.hasSettledThisHour()) {
                            // Save arrival
                            if (rp.recordArrivals) {
                                IOUtils.arrivalToFile(part, currentIsoDate, currentHour, "arrivals_" + today + ".dat", true);
                            }
                            // Add arrival to connectivity file
                            int destIndex = siteEndNames.indexOf(part.getLastArrival()); // TODO: make sure arrival indexes are for siteEnd, not site
                            connectivity[part.getStartIndex()][destIndex] += part.getDensity();

                            // Reset ability to settle
                            part.setSettledThisHour(false);
                        }
                    }

                    // Hourly updates to pstep arrays
                    if (rp.recordPsteps) {
                        IOUtils.pstepsUpdater(particles, rp, pstepsMature, pstepsImmature, 3600);
                    }

                    // Write Psteps
                    if (rp.recordPsteps && stepcount % (rp.pstepsInterval * rp.stepsPerStep) == 0) {
                        // Trim arrays to non-zero rows and write to file
                        float[][] psImmTrim = null;
                        try {
                            psImmTrim = nonZeroRows(pstepsImmature);
                        } catch (Exception ignored) {
                        }
                        float[][] psMatTrim = null;
                        try {
                            psMatTrim = nonZeroRows(pstepsMature);
                        } catch (Exception ignored) {
                        }
                        System.out.println("Writing psteps");
                        if (psImmTrim != null) {
                            IOUtils.writeFloatArrayToFile(psImmTrim, "pstepsImmature_" + today + "_" + stepcount + ".dat", false, true);
                        }
                        if (psMatTrim != null) {
                            IOUtils.writeFloatArrayToFile(psMatTrim, "pstepsMature_" + today + "_" + stepcount + ".dat", false, true);
                        }

                        pstepsImmature = new float[meshes.get(0).getNElems()][habitat.size()];
                        pstepsMature = new float[meshes.get(0).getNElems()][habitat.size()];
                    }

                    if (rp.recordConnectivity && stepcount % (rp.connectivityInterval * rp.stepsPerStep) == 0) {
                        IOUtils.writeFloatArrayToFile(connectivity, "connectivity_" + today + "_" + stepcount + ".dat", false, false);
                        connectivity = new float[habitat.size()][habitatEnd.size()];
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
            IOUtils.printFileHeader(particleRestartHeader, "locationsEnd_" + currentIsoDate.getDateStr() + ".dat");
            IOUtils.particlesToRestartFile(particles, 0, "locationsEnd_" + currentIsoDate.getDateStr() + ".dat", true, rp, rp.nparts * rp.numberOfDays * 10);
            if (rp.recordElemActivity) {
                IOUtils.writeIntegerArrayToFile(elemActivity, "elementActivity.dat");
            }

            executorService.shutdownNow();
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            executorService.shutdownNow();
        }

        FileWriter fstream = new FileWriter("startSitesUsed.dat", false);
        PrintWriter out = new PrintWriter(fstream);
        out.println("site\tx\ty\tdepth\tmesh\tcentroid\telem\tMeshType\tu\tv\tw\tsalinity");
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
                                                    ISO_datestr currentDate, int currentTime, double startDensity, int numParticlesCreated) {

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

            Particle p = new Particle(xstart, ystart, rp.startDepth, habitat.get(startid).getID(), startid, numParticlesCreated + i,
                    rp.mortalityRate, startDensity, currentDate, currentTime, rp.coordRef, rp.species);
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
    public static List<HydroField> readHydroFields(List<Mesh> meshes, ISO_datestr currentIsoDate, int currentHour, boolean isLastDay, RunProperties rp) {
        List<HydroField> hydroFields = new ArrayList<>();

        // 24 hr files only case - read once a day
        int m = 0;
        for (Mesh mesh : meshes) {
            if (mesh.getType().equalsIgnoreCase("FVCOM")) {

                try {
                    System.out.println("Reading file " + currentHour); // Dima file naming format: minch2_20171229_0003.nc
                    String[] varNames1 = new String[]{"u", "v", "ww", "salinity", "temp", "zeta", "km", "short_wave"};

                    // Normal "forwards time"
                    if (!rp.backwards) {
                        ISO_datestr tomorrow = ISO_datestr.getTomorrow(currentIsoDate);
                        if (isLastDay) {
                            System.out.println("** Last day - reading same hydro file twice **");
                            tomorrow = currentIsoDate;
                        }
                        String[] dirsT1 = new String[]{
                                rp.datadir + rp.datadirPrefix + currentIsoDate.getYear() + rp.datadirSuffix + System.getProperty("file.separator"),
                                rp.datadir2 + rp.datadir2Prefix + currentIsoDate.getYear() + rp.datadir2Suffix + System.getProperty("file.separator")};
                        String[] dirsT2 = new String[]{
                                rp.datadir + rp.datadirPrefix + tomorrow.getYear() + rp.datadirSuffix + System.getProperty("file.separator"),
                                rp.datadir2 + rp.datadir2Prefix + tomorrow.getYear() + rp.datadir2Suffix + System.getProperty("file.separator")};
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
                                new File(rp.datadir + rp.datadirPrefix + currentIsoDate.getYear() + rp.datadirSuffix + System.getProperty("file.separator")),
                                new WildcardFileFilter(rp.location + rp.minchVersion + "_" + currentIsoDate.getYear() + String.format("%02d", currentIsoDate.getMonth()) + String.format("%02d", currentIsoDate.getDay()) + "*_rev.nc"),
                                null);
                        ISO_datestr yesterday = ISO_datestr.getYesterday(currentIsoDate);
                        if (isLastDay) {
                            System.out.println("** Last day - reading same hydro file twice **");
                            yesterday = currentIsoDate;
                        }
                        List<File> files2 = (List<File>) FileUtils.listFiles(
                                new File(rp.datadir + rp.datadirPrefix + yesterday.getYear() + rp.datadirSuffix + System.getProperty("file.separator")),
                                new WildcardFileFilter(rp.location + rp.minchVersion + "_" + yesterday.getYear() + String.format("%02d", yesterday.getMonth()) + String.format("%02d", yesterday.getDay()) + "*_rev.nc"),
                                null);
                        // Read both files and combine
                        hydroFields.add(new HydroField(files1.get(0).getCanonicalPath(), files2.get(0).getCanonicalPath(), varNames1, null, null, null, "FVCOM", rp));
                    }

                } catch (Exception e) {
                    System.out.println("Hydro file not found, check PROPERTIES: datadir, datadirPrefix, datadirSuffix, location, minchVersion");
                    if (!rp.backwards) {
                        System.err.println("Requested file: " + rp.datadir + rp.datadirPrefix + currentIsoDate.getYear() + rp.datadirSuffix + System.getProperty("file.separator")
                                + rp.location + rp.minchVersion + "_" + currentIsoDate.getYear() + String.format("%02d", currentIsoDate.getMonth()) + String.format("%02d", currentIsoDate.getDay()) + "*.nc");
                    } else {
                        System.err.println("Requested file: " + rp.datadir + rp.datadirPrefix + currentIsoDate.getYear() + rp.datadirSuffix + System.getProperty("file.separator")
                                + rp.location + rp.minchVersion + "_" + currentIsoDate.getYear() + String.format("%02d", currentIsoDate.getMonth()) + String.format("%02d", currentIsoDate.getDay()) + "*_rev.nc");
                    }
                    System.exit(1);
                }
            } else if (mesh.getType().equalsIgnoreCase("ROMS_TRI")) {
                String filename1 = rp.datadir2 + rp.datadir2Prefix + currentIsoDate.getYear() + rp.datadir2Suffix + System.getProperty("file.separator")
                        + "NEATL_" + currentIsoDate.getYear() + String.format("%02d", currentIsoDate.getMonth()) + String.format("%02d", currentIsoDate.getDay()) + ".nc";
                String filename2 = rp.datadir2 + rp.datadir2Prefix + currentIsoDate.getYear() + rp.datadir2Suffix + System.getProperty("file.separator")
                        + "NEATL_" + currentIsoDate.getYear() + String.format("%02d", currentIsoDate.getMonth()) + String.format("%02d", currentIsoDate.getDay()) + ".nc";
                String[] varNames1 = {"u", "v", "", "", ""};
                // Read both files and combine
                hydroFields.add(new HydroField(filename1, filename2, varNames1, null, null, null, "ROMS_TRI", rp));
            }

            m++;

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
        System.out.println("count " + count);
        float[][] temp = null;
        if (count > 0) {
            temp = new float[count][A[0].length + 1];
            System.out.println("temp size = " + temp.length + " " + temp[0].length);
            //System.out.println("A size = "+A.length+" "+A[0].length);
            int p = 0;
            for (int row : list) {
                temp[p][0] = row;
                for (int j = 0; j < A[0].length; j++) {
                    //System.out.println(A[row][j]);
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
