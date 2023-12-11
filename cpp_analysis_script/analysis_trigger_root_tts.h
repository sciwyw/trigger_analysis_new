//
// Created by ineffablord on 23-11-28.
//

#ifndef ANALYSIS_K40_ANALYSIS_TRIGGER_ROOT_TTS_H
#define ANALYSIS_K40_ANALYSIS_TRIGGER_ROOT_TTS_H

struct QEdata{
    Double_t wl;
    Double_t qe;
};
Double_t interpolate(const std::vector<QEdata>& , Double_t );

#endif //ANALYSIS_K40_ANALYSIS_TRIGGER_ROOT_TTS_H
