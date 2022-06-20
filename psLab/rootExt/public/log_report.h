#ifndef ROOTEXT_LOG_REPORT_H_
#define ROOTEXT_LOG_REPORT_H_

//#define log_fatal(format, ...) {printf(format, ##__VA_ARGS__); StopOrGo();}
// exit(1) is better solution if we are running in batch mode
#define log_fatal(format, ...) {printf(format, ##__VA_ARGS__); exit(1);}
#define log_error(format, ...) printf(format, ##__VA_ARGS__)
#define log_warn(format, ...) printf(format, ##__VA_ARGS__)
#define log_info(format, ...) printf(format, ##__VA_ARGS__)
#define log_debug(format, ...) printf(format, ##__VA_ARGS__)

#endif // ROOTEXT_LOG_REPORT_H_
