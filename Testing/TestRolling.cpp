#include <iostream>
#include <map>

#include "RollingSingles.h"
#include "ReadCSV.h"

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <csv_or_bin_file>" << std::endl;
        return 1;
    }

    const std::string filename = argv[1];
    double duration = 0.0;
    auto chunk = readFileAuto(filename, duration);

    RollingSingles rolling(400);
    rolling.appendChunk(chunk);
    std::cout << "Total measurement time " << duration << std::endl;

    for (const auto &entry : rolling.allChannels()) {
        const Singles &s = entry.second;
        std::cout << "Channel " << entry.first << " baseSecond=" << s.baseSecond
                  << " buckets=" << s.eventsPerSecond.size() << std::endl;
    }

    return 0;
}
