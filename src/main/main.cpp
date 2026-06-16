#include "sim/simulation_config.h"
#include "sim/simulation_engine.h"
#include "sim/json_config.h"
#include "sysbin/dmdin.h"
#include <iostream>
#include <string>
#include <stdexcept>
#include <cstring>

static void print_usage(const char* prog) {
    std::cerr << "Usage:\n"
              << "  " << prog << " run <system.dmdin> <config.json>\n"
              << "  " << prog << " config --template\n"
              << "  " << prog << " --help\n";
}

static int cmd_run(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Error: missing arguments.\n";
        print_usage(argv[0]);
        return 1;
    }
    std::string dmdin_path = argv[2];
    std::string config_path = argv[3];

    SimulationConfig cfg = read_dmdin(dmdin_path);
    apply_json_config(cfg, config_path);

    SimulationEngine engine = build_simulation(cfg);
    engine.run();

    return 0;
}

static int cmd_config(int argc, char* argv[]) {
    if (argc < 3 || std::strcmp(argv[2], "--template") != 0) {
        std::cerr << "Error: unknown option.\n";
        print_usage(argv[0]);
        return 1;
    }
    std::cout << generate_config_template() << "\n";
    return 0;
}

int main(int argc, char* argv[]) {
    try {
        if (argc < 2) {
            print_usage(argv[0]);
            return 1;
        }
        if (std::strcmp(argv[1], "--help") == 0) {
            print_usage(argv[0]);
            return 0;
        }
        if (std::strcmp(argv[1], "run") == 0) {
            return cmd_run(argc, argv);
        }
        if (std::strcmp(argv[1], "config") == 0) {
            return cmd_config(argc, argv);
        }
        std::cerr << "Error: unknown command '" << argv[1] << "'.\n";
        print_usage(argv[0]);
        return 1;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
