/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 *
 * This file is based on a class given on
 * https://stackoverflow.com/questions/4010185/parallel-for-for-java
 *
 */
package particle_track;

import java.util.Collection;
import java.util.List;
import java.util.ArrayList;
import java.util.concurrent.Callable;
import java.util.concurrent.Executors;
import java.util.concurrent.ExecutorService;
import java.lang.InterruptedException;

/**
 *
 * @author sa01ta
 */
public class Parallel {
    private static final int NUM_CORES = Runtime.getRuntime().availableProcessors();

    private static final ExecutorService forPool = Executors.newFixedThreadPool(NUM_CORES * 2, Executors.defaultThreadFactory());

    public static <T> void For(final Iterable<T> elements, final Operation<T> operation) {
        try {
            // invokeAll blocks for us until all submitted tasks in the call complete
            forPool.invokeAll(createCallables(elements, operation));
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }

    public static <T> Collection<Callable<Void>> createCallables(final Iterable<T> elements, final Operation<T> operation) {
        List<Callable<Void>> callables = new ArrayList<Callable<Void>>();
        for (final T elem : elements) {
            callables.add(new Callable<Void>() {
                @Override
                public Void call() {
                    operation.perform(elem);
                    return null;
                }
            });
        }

        return callables;
    }

    public static interface Operation<T> {
        public void perform(T pParameter);
    }
}

//// How do we make the above actually do something?
//// Collection of items to process in parallel
//Collection<Integer> elems = new ArrayList<Integer>();
//for (int i = 0; i < 40; ++i) {
//    elems.add(i);
//}
//Parallel.For(elems, 
// // The operation to perform with each item
// new Parallel.Operation<Integer>() {
//    public void perform(Integer param) {
//        System.out.println(param);
//    };
//});


// Basic process given on
// https://stackoverflow.com/questions/4010185/parallel-for-for-java
//    ExecutorService exec = Executors.newFixedThreadPool(SOME_NUM_OF_THREADS);
//try {
//    for (final Object o : list) {
//        exec.submit(new Runnable() {
//            @Override
//            public void run() {
//                // do stuff with o.
//            }
//        });
//    }
//} finally {
//    exec.shutdown();
//}

    

