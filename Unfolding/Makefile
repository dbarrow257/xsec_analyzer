all: Unfolder Validation

Unfolder: Unfolder.C 
	$(CXX) $(shell root-config --cflags --libs) -O3 -I$(STV_ANALYSIS_DIR)/Utils/Includes $(STV_ANALYSIS_DIR)/Bin/*.o -o $@ $^

Validation: Validation.C 
	$(CXX) $(shell root-config --cflags --libs) -O3 -I$(STV_ANALYSIS_DIR)/Utils/Includes $(STV_ANALYSIS_DIR)/Bin/*.o -o $@ $^

clean:
	$(RM) Unfolder Validation
	$(RM) *~
	$(RM) stv_root_dict_rdict.pcm
