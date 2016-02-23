#ifndef H_ADC_PV
#define H_ADC_PV
typedef const unsigned int addr_t;

#define ADC_OUTPUT_PREC 12

addr_t ADC_BASE		= 0x40038000;
addr_t ADC_CR		= ADC_BASE + 0;
addr_t ADC_MR		= ADC_BASE + 0x04;
addr_t ADC_CHER		= ADC_BASE + 0x08;
addr_t ADC_CHDR		= ADC_BASE + 0x10;
addr_t ADC_CHSR		= ADC_BASE + 0x18;
addr_t ADC_CGR		= ADC_BASE + 0x48;
addr_t ADC_COR		= ADC_BASE + 0x4C;
addr_t ADC_CDR0		= ADC_BASE + 0x50;
addr_t ADC_CDR1		= ADC_BASE + 0x54;
addr_t ADC_CDR2		= ADC_BASE + 0x58;
addr_t ADC_CDR3		= ADC_BASE + 0x5C;
addr_t ADC_CDR4		= ADC_BASE + 0x60;
addr_t ADC_CDR5		= ADC_BASE + 0x64;
addr_t ADC_CDR6		= ADC_BASE + 0x68;
addr_t ADC_CDR7		= ADC_BASE + 0x6C;
addr_t ADC_CDR8		= ADC_BASE + 0x70;
addr_t ADC_CDR9		= ADC_BASE + 0x74;
addr_t ADC_CDR10	= ADC_BASE + 0x78;
addr_t ADC_CDR11	= ADC_BASE + 0x7C;
addr_t ADC_CDR12	= ADC_BASE + 0x80;
addr_t ADC_CDR13	= ADC_BASE + 0x84;
addr_t ADC_CDR14	= ADC_BASE + 0x88;
addr_t ADC_CDR15	= ADC_BASE + 0x8C;
addr_t ADC_WPMR		= ADC_BASE + 0xE4;

const unsigned int ADC_PAYLOAD_LENGTH = 4;


extern void ADC_Configuration(unsigned int addr, int value);

extern unsigned int ADC_GetStatus();

extern unsigned int ADC_GetData(unsigned short index);





#endif