
#pragma warning(disable : 4996)
#ifndef FILE_CONTROL_PARTICLE
#define FILE_CONTROL_PARTICLE
#include "D:\VisualStudioProj\My_Class\Particle_Music\Particle.h"

//文件格式：每一行分别为：起始时间	（水平制表符） 音高


FILE* PARTICLEFILECONTROL = nullptr;
FILE* SBPARTICLEFILECONTROL = nullptr;
char PARTICLEFILELOCATE[255] = "Input.txt";
char SBPARTICLEFILELOCATE[255] = "Input_SB.txt";
vector<vector<int>> MIDI_INFO_PARTICLE;

void open_file() {
	PARTICLEFILECONTROL = fopen(PARTICLEFILELOCATE, "a+");
}

void load_txt_file(){
	int TEMP_NOTE[2];
	for (int i=0; fscanf(PARTICLEFILECONTROL, "%d %d", TEMP_NOTE, TEMP_NOTE + 1)!=EOF;) {
		if (MIDI_INFO_PARTICLE.size() == 0) {
			vector<int> TEM_NOTE;
			TEM_NOTE.push_back(*TEMP_NOTE);
			TEM_NOTE.push_back(*(TEMP_NOTE + 1));
			MIDI_INFO_PARTICLE.push_back(TEM_NOTE);
		}
		else {
			if (*TEMP_NOTE == MIDI_INFO_PARTICLE[i][0]) {
				MIDI_INFO_PARTICLE[i].push_back(*(TEMP_NOTE + 1));
			}
			else {
				i++;
				vector<int> TEM_NOTE;
				TEM_NOTE.push_back(*TEMP_NOTE);
				TEM_NOTE.push_back(*(TEMP_NOTE + 1));
				MIDI_INFO_PARTICLE.push_back(TEM_NOTE);
			}
		}
	}
}

void print_txt_file() {
	for (int i = 0; i < (int)MIDI_INFO_PARTICLE.size(); i++) {
		for (int j = 0; j < (int)MIDI_INFO_PARTICLE[i].size(); j++) {
			printf("%d    ", MIDI_INFO_PARTICLE[i][j]);
		}
		printf("\n");
	}
}

void open_sound_block_file() {
	SBPARTICLEFILECONTROL = fopen(SBPARTICLEFILELOCATE, "a+");
}

vector <vector<int>> SBNOTEINFO;

void load_SB_Particle_file() {
	srand(time(NULL));
	int SBTEMP_NOTE[3];
	for (int i = 0; fscanf(SBPARTICLEFILECONTROL, "%d %d %d", SBTEMP_NOTE, SBTEMP_NOTE + 1, SBTEMP_NOTE + 2) != EOF;) {
		vector<int> TEMPSB;
		TEMPSB.push_back(*SBTEMP_NOTE);
		TEMPSB.push_back(*(SBTEMP_NOTE+1));
		TEMPSB.push_back(*(SBTEMP_NOTE+2));
		SBNOTEINFO.push_back(TEMPSB);
	}
}

void print_SB_file() {
	for (int i = 0; i < (int)SBNOTEINFO.size(); i++) {
		for (int j = 0; j < 3; j++) {
			printf("%d    ", SBNOTEINFO[i][j]);
		}
		printf("\n");
	}
}

#endif