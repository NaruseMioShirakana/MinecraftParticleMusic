#ifndef MCPARTICLEMIDIPLAYERMAINCONTROL_
#define MCPARTICLEMIDIPLAYERMAINCONTROL_
#pragma warning(disable : 4996)
#include "D:\VisualStudioProj\My_Class\Particle_Music\Particle_file_control.h"
int StartPos[2];
int TempYPos = 106;
char ParMode = 'a';
int Tem_SCANF;
//-50 51
char TemArr_SB[255] = "minecraft:sea_latern\0";
void Tran_MODEA() {//A
	printf("\n\n\n\n������һ�е�������꣨x z����\n");
	Tem_SCANF = Tem_SCANF = scanf("%d %d", StartPos, StartPos+1);

	for (int i = 0; i < (int)MIDI_INFO_PARTICLE.size(); i++) {
		for (int j = 1; j < (int)MIDI_INFO_PARTICLE[i].size(); j++) {
			if (i == 0) {
				Butter_Fly_To_Next_Point_Sim(*StartPos, MIDI_INFO_PARTICLE[i][0], TempYPos, TempYPos, *(StartPos + 1), MIDI_INFO_PARTICLE[i][j], 20, "Tem");//���˺���������һ��else�к�������ΪParticle.h�е��κ�һ���������������з���������Բ�ĺ����ɸı���������
			}
			else {
				Butter_Fly_To_Next_Point_Sim(MIDI_INFO_PARTICLE[i-1][0], MIDI_INFO_PARTICLE[i][0], TempYPos, TempYPos, MIDI_INFO_PARTICLE[i-1][1], MIDI_INFO_PARTICLE[i][j], 20, "Tem");
			}
		}
	}
	printf("ת������\n");
	return;
}

void Tran_MODEB() {//B
	printf("\n\n\n\n������һ�е�������꣨x z����\n");
	Tem_SCANF = Tem_SCANF = scanf("%d %d", StartPos, StartPos + 1);

	for (int i = 0; i < (int)MIDI_INFO_PARTICLE.size(); i++) {
		for (int j = 1; j < (int)MIDI_INFO_PARTICLE[i].size(); j++) {
			if (i == 0) {
				straight_line(*StartPos, MIDI_INFO_PARTICLE[i][0], TempYPos, TempYPos, *(StartPos + 1), MIDI_INFO_PARTICLE[i][j], 20, "Tem");//���˺���������һ��else�к�������ΪParticle.h�е��κ�һ���������������з���������Բ�ĺ����ɸı���������
			}
			else {
				if (i > (int)MIDI_INFO_PARTICLE[i - 1].size()) {
					straight_line(MIDI_INFO_PARTICLE[i - 1][0], MIDI_INFO_PARTICLE[i][0], TempYPos, TempYPos, MIDI_INFO_PARTICLE[i - 1][(int)MIDI_INFO_PARTICLE[i - 1].size()-1], MIDI_INFO_PARTICLE[i][j], 20, "Tem");
				}
				else {
					straight_line(MIDI_INFO_PARTICLE[i - 1][0], MIDI_INFO_PARTICLE[i][0], TempYPos, TempYPos, MIDI_INFO_PARTICLE[i - 1][j], MIDI_INFO_PARTICLE[i][j], 20, "Tem");
				}
				
			}
		}
	}
	printf("ת������\n");
	return;
}

void Tran_MODEC() {//C
	printf("\n\n\n\n������һ�е�������꣨x z����\n");
	Tem_SCANF = Tem_SCANF = scanf("%d %d", StartPos, StartPos + 1);

	for (int i = 0; i < (int)MIDI_INFO_PARTICLE.size(); i++) {
		for (int j = 1; j < (int)MIDI_INFO_PARTICLE[i].size(); j++) {
			if (i == 0) {
				straight_line(*StartPos, MIDI_INFO_PARTICLE[i][0], TempYPos, TempYPos, *(StartPos + 1), MIDI_INFO_PARTICLE[i][j], 20, "Tem");//���˺���������һ��else�к�������ΪParticle.h�е��κ�һ���������������з���������Բ�ĺ����ɸı���������
			}
			else {
				straight_line(MIDI_INFO_PARTICLE[i - 1][0], MIDI_INFO_PARTICLE[i][0], TempYPos, TempYPos, MIDI_INFO_PARTICLE[i - 1][1], MIDI_INFO_PARTICLE[i][j], 20, "Tem");
			}
		}
		if ((int)MIDI_INFO_PARTICLE[i - 1].size() > (int)MIDI_INFO_PARTICLE[i].size()) {
			for (int j = (int)MIDI_INFO_PARTICLE[i].size() - 1; j < (int)MIDI_INFO_PARTICLE[i - 1].size(); j++) {
				straight_line(MIDI_INFO_PARTICLE[i - 1][0], MIDI_INFO_PARTICLE[i][0], TempYPos, TempYPos, MIDI_INFO_PARTICLE[i - 1][j], MIDI_INFO_PARTICLE[i][1], 20, "Tem");
			}
		}
	}
	printf("ת������\n");
	return;
}

void Tran_MODED() {//D
	printf("\n\n\n\n������һ�е�������꣨x z����\n");
	Tem_SCANF = Tem_SCANF = scanf("%d %d", StartPos, StartPos + 1);

	for (int i = 0; i < (int)MIDI_INFO_PARTICLE.size(); i++) {
		for (int j = 1; j < (int)MIDI_INFO_PARTICLE[i].size(); j++) {
			if (i == 0) {
				straight_line(*StartPos, MIDI_INFO_PARTICLE[i][0], TempYPos, TempYPos, *(StartPos + 1), MIDI_INFO_PARTICLE[i][j], 20, "Tem");//���˺���������һ��else�к�������ΪParticle.h�е��κ�һ���������������з���������Բ�ĺ����ɸı���������
			}
			else {
				if (i > (int)MIDI_INFO_PARTICLE[i - 1].size()) {
					straight_line(MIDI_INFO_PARTICLE[i - 1][0], MIDI_INFO_PARTICLE[i][0], TempYPos, TempYPos, MIDI_INFO_PARTICLE[i - 1][(int)MIDI_INFO_PARTICLE[i - 1].size() - 1], MIDI_INFO_PARTICLE[i][j], 20, "Tem");
				}
				else {
					straight_line(MIDI_INFO_PARTICLE[i - 1][0], MIDI_INFO_PARTICLE[i][0], TempYPos, TempYPos, MIDI_INFO_PARTICLE[i - 1][j], MIDI_INFO_PARTICLE[i][j], 20, "Tem");
				}
			}
		}
		if ((int)MIDI_INFO_PARTICLE[i - 1].size() > (int)MIDI_INFO_PARTICLE[i].size()) {
			for (int j = (int)MIDI_INFO_PARTICLE[i].size() - 1; j < (int)MIDI_INFO_PARTICLE[i - 1].size(); j++) {
				straight_line(MIDI_INFO_PARTICLE[i - 1][0], MIDI_INFO_PARTICLE[i][0], TempYPos, TempYPos, MIDI_INFO_PARTICLE[i - 1][j], MIDI_INFO_PARTICLE[i][(int)MIDI_INFO_PARTICLE[i].size() - 1], 20, "Tem");
			}
		}
	}
	printf("ת������\n");
	return;
}
void StartProject() {
	file1 = fopen("Temp.mcfunction", "a+");//�޸��ļ���
	open_file();
	load_txt_file();
	printf("ÿһ�м�ΪͬһT��������\n");
	print_txt_file();
	PARTICLE_ERROR0_:
	printf("\n\n\n\n\n������y���꣨>0)��y = ");
	Tem_SCANF = scanf("%d", &TempYPos);
	if (TempYPos <= 0) {
		printf("�����ʽ����");
		goto PARTICLE_ERROR0_;
	}
	else {
		PARTICLE_ERROR1_:
		printf("\n\n\n\n\n������ģʽ��ABCD����\nģʽA��ÿһ�η�֧�����������ۺϣ���ÿһ�η�֧�����Ϊ����ĵ�һ�����������꣩\nģʽB��ÿһ�η�֧�����������ۺϣ���ÿһ�η�֧�����Ϊ���������ͬ��������������������Ŀ������һ�У��������ִӶ�������һ����ʼ��\nģʽC��ÿһ�η�֧����������ۺϣ��ۺϵ����ĵ�һ�����������꣩\nģʽD��ÿһ�η�֧����������ۺϣ��ۺϵ����������Bģʽ��\nģʽ = ");
		Tem_SCANF = scanf("%c", &ParMode);//���ڻ�������������һ������ĵĻس�����
		Tem_SCANF = scanf("%c", &ParMode);
		switch (ParMode) {
		case 'A': Tran_MODEA(); break;
		case 'B': Tran_MODEB(); break;
		case 'C': Tran_MODEC(); break;
		case 'D': Tran_MODED(); break;
		default: printf("�����������дA B C D�е�һ��\n"); goto PARTICLE_ERROR1_; break;
		}
	}
	fclose(file1);
	return;
}

void Tran_MODEROUND() {//ROUND
	for (int i = 1; i < (int)MIDI_INFO_PARTICLE.size(); i++) {
		if ((int)MIDI_INFO_PARTICLE[i].size() > (int)MIDI_INFO_PARTICLE[i - 1].size()) {
			Init_Reset((int)MIDI_INFO_PARTICLE[i - 1].size());
		}
		for (int j = 1; j < (int)MIDI_INFO_PARTICLE[i].size(); j++) {
			if (j < (int)MIDI_INFO_PARTICLE[i - 1].size()) {
				Round_Truck_lisan(MIDI_INFO_PARTICLE[i - 1][0], MIDI_INFO_PARTICLE[i][0], TempYPos, TempYPos, MIDI_INFO_PARTICLE[i - 1][j], MIDI_INFO_PARTICLE[i][j], j);
			}
			else {
				Round_Truck_lisan(MIDI_INFO_PARTICLE[i - 1][0], MIDI_INFO_PARTICLE[i][0], TempYPos, TempYPos, MIDI_INFO_PARTICLE[i - 1][(int)MIDI_INFO_PARTICLE[i - 1].size()-1], MIDI_INFO_PARTICLE[i][j], j);
			}
		}
	}
	printf("ת������\n");
	return;
}

void StartRound() {
	file1 = fopen("Temp.mcfunction", "a+");//�޸��ļ���
	open_file();
	load_txt_file();
	printf("ÿһ�м�ΪͬһT��������\n");
	print_txt_file();
	PARTICLE_ERROR2_:
	printf("\n\n\n\n\n������y���꣨>0)��y = ");
	Tem_SCANF = scanf("%d", &TempYPos);
	if (TempYPos <= 0) {
		printf("�����ʽ����");
		goto PARTICLE_ERROR2_;
	}
	else {
		Tran_MODEROUND();
	}
	fclose(file1);
	return;
}

//execute if score @p Timer matches 5385 run execute as @p at @p run playsound flp.80 player @a ~ ~ ~ 0.64 1

void playsound_setblock() {
	file1 = fopen("Temp.mcfunction", "a+");
	open_sound_block_file();
	load_SB_Particle_file();
	print_SB_file();
	for (int i = 0; i < (int)SBNOTEINFO.size(); i++) {
		fprintf(file1, "execute if score @p Timer matches %d run execute as @p at @p run playsound flp.%d player @a ~ ~ ~ %.2f 1\n", SBNOTEINFO[i][0], SBNOTEINFO[i][1], (double)SBNOTEINFO[i][2] / 100.0);
		fprintf(file1, "execute if score @p Timer matches %d run setblock %d %d %d %s replace\n", SBNOTEINFO[i][0], SBNOTEINFO[i][0], TempYPos, SBNOTEINFO[i][1], TemArr_SB);
	}
	fclose(file1);
}

void EndParticle() {
	file1 = fopen("Temp.mcfunction", "a+");
	open_sound_block_file();
	load_SB_Particle_file();
	print_SB_file();
	for (int i = 0; i < (int)SBNOTEINFO.size(); i++) {
		Butter_Fly_Spread(SBNOTEINFO[i][0], TempYPos, SBNOTEINFO[i][1],1,1.0);
	}
	fclose(file1);
}



























#endif // !MCPARTICLEMIDIPLAYERMAINCONTROL_
