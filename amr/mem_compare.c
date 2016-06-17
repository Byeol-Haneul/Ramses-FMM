#include <stdint.h>
uint32_t  memcmp_(void *str1, void *str2, uint32_t *n)
{
  return memcmp(str1, str2, *n);
}
