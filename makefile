simple_integration: simple_integration.cc ./include* 
	/opt/homebrew/bin/h5c++ -o simple_integration simple_integration.cc -I./include

simple_integration_frames: simple_integration_frames.cc ./include* 
	/opt/homebrew/bin/h5c++ -o simple_integration_frames simple_integration_frames.cc -I./include
