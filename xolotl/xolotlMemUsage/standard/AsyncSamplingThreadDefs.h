#ifndef XMEMUSAGE_ASYNC_SAMPLING_THREAD_DEFS_H
#define XMEMUSAGE_ASYNC_SAMPLING_THREAD_DEFS_H

// This file contains implementations of AsyncSamplingThread template 
// class member functions.  It is needed to break a dependency cycle
// between AsyncSampler and AsyncSamplingThread.
#include "AsyncSamplingThread.h"

namespace xolotlMemUsage {

template<typename RPD, typename PSD>
void
AsyncSamplingThread<RPD, PSD>::StartSamplingFor(AsyncSampler<RPD, PSD>* sampler)
{
    // Ensure we are the only ones updating our data structures.
    std::lock_guard<std::mutex> ourLock(mtx);

    // Ensure we have a sampling thread.
    if(not haveSamplingThread)
    {
        // We do not have a sampling thread.  Create one.
        haveSamplingThread = true;
        samplingThreadDone = false;
        samplingThreadResult = std::async(&AsyncSamplingThread::CollectSamples, this);
    }

    // Ensure we know about the sampler
    auto arIter = allSamplers.find(sampler->GetName());
    if(arIter == allSamplers.end())
    {
        // We did not already know about a sampler with this name.
        allSamplers.emplace(sampler->GetName(), sampler);
    }
    else if(arIter->second != sampler)
    {
        // We knew about a sampler with this name, but it isn't 
        // the one that we have been given.
        std::cerr << "WARNING: unmatched start/stop sampling for async samplers."
            << "Replacing sampler with name " << sampler->GetName() 
            << ".  Old sampler's sampled data will be lost." 
            << std::endl;
        auto oldSampler = arIter->second;
        allSamplers[sampler->GetName()] = sampler;

        // Ensure that the old sampler is no longer considered active.
        auto activeIter = activeSamplers.find(oldSampler);
        if(activeIter != activeSamplers.end())
        {
            // Old sampler object was active, so remove it so that new
            // sampler object can safely take its place.
            activeSamplers.erase(activeIter);
        }
    }

    // Ensure the sampler is considered active.
    activeSamplers.emplace(sampler);
}


template<typename RPD, typename PSD>
void
AsyncSamplingThread<RPD, PSD>::StopSamplingFor(AsyncSampler<RPD, PSD>* sampler)
{
    // Ensure we are the only ones updating our data structures.
    std::lock_guard<std::mutex> ourLock(mtx);

    // Remove the sampler from the set of active samplers.
    auto activeIter = activeSamplers.find(sampler);
    if(activeIter != activeSamplers.end())
    {
        activeSamplers.erase(activeIter);

        // Check whether we have *any* samplers that are active.
        if(activeSamplers.empty())
        {
            // We have no active samplers - no need
            // to keep the sampling thread running.
            samplingThreadDone = true;        
            samplingThreadResult.get();   // wait till sample thread goes away
            haveSamplingThread = false;
        }
    }
}

template<typename RPD, typename PSD>
void
AsyncSamplingThread<RPD, PSD>::CollectSamples(void)
{
    while(not samplingThreadDone)
    {
        Sample();
        std::this_thread::sleep_for(samplingInterval);
    }
}

} // end namespace xolotlMemUsage

#endif // XMEMUSAGE_ASYNC_SAMPLING_THREAD_DEFS_H
