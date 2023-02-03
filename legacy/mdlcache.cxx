#include "mdl.h"
#include <array>
#include <memory>
typedef std::array<int64_t,NDIM+1> KEY;
typedef mdl::hash::HASH<KEY> HASH;
typedef mdl::ARC<KEY> ARC;

// Create a hash table that will hold at most "nMaxElements" entries
// Each key will be KeyLength long in units of "uint32_t"
extern "C"
void *HashCreate(int nMaxElements) {
    return new HASH(nMaxElements);
    }

extern "C"
void HashDestroy(void *pTable) {
    delete reinterpret_cast<HASH*>(pTable);
    }

extern "C"
void HashClear(void *pTable) {
    reinterpret_cast<HASH*>(pTable)->clear();
    }

extern "C"
void *HashLookup(void *pTable, uint32_t uHash,const KEY &key) {
    return reinterpret_cast<HASH*>(pTable)->lookup(uHash,key);
    }

extern "C"
void HashInsert(void *pTable, uint32_t uHash,const KEY &key,void *data) {
    return reinterpret_cast<HASH*>(pTable)->insert(uHash,key,data);
    }

extern "C"
void HashRemove(void *pTable, uint32_t uHash,const KEY &key) {
    reinterpret_cast<HASH*>(pTable)->remove(uHash,key);
    }

extern "C"
void HashStatistics(void *pTable) {
    reinterpret_cast<HASH*>(pTable)->print_statistics();
    }

using namespace mdl;
class ramsesCACHEhelper : public CACHEhelper {
public:
    typedef void (*pack_function_type)   (const void *src,const uint32_t &size,      void *dst);
    typedef void (*unpack_function_type) (void       *dst,const uint32_t &size,const void *src, const void *key);
    typedef void (*init_function_type)   (void       *dst);
protected:
    uint32_t nPack, nFlush;
    void *ctx;
    pack_function_type   pack_function;
    unpack_function_type unpack_function;
    init_function_type   init_function;
    pack_function_type   flush_function;
    unpack_function_type combine_function;
    uint32_t (*get_thread)  (void *,uint32_t,uint32_t,const void *);
    void* (*create_function) (void *,uint32_t,const void *);
    uint32_t (*hash_func)(const void *key);
    int32_t (*get_tile)(void *ctx,const void *data,uint32_t nkey, void * const *keys,void const **grid);

    virtual void    pack(void *dst, const void *src)                  override { (*pack_function)(src,nPack/4,dst); }
    virtual void  unpack(void *dst, const void *src, const void *key) override { (*unpack_function)(dst,nPack/4,src,key); }
    virtual void    init(void *dst)                                   override { if (init_function)   (*init_function)(dst); }
    virtual void   flush(void *dst, const void *src)                  override { if (flush_function)  (*flush_function)(src,nPack/4,dst); }
    virtual void combine(void *dst, const void *src, const void *key) override { if (combine_function)(*combine_function)(dst,nPack/4,src,key); }
    virtual void *create(uint32_t size, const void *pKey) override
	{ return init_function==nullptr ? (*create_function)(ctx,size,pKey) : nullptr; }
    virtual uint32_t getThread(uint32_t uLine, uint32_t uId, uint32_t size, const void *pKey) override
	{ return (*get_thread)(ctx,uLine,size,pKey); }
    virtual uint32_t getTile(void const * data, uint32_t key_size, void const *pKey, uint32_t n, void * const *keys, void const **tile) override
	{     return (*get_tile)(ctx,data,n,keys,tile); }
    virtual uint32_t hash(const void *key) override {return (*hash_func)(key);}
    virtual uint32_t pack_size()  override {return nPack;}
    virtual uint32_t flush_size() override {return nFlush;}
public:
    explicit ramsesCACHEhelper(uint32_t nData,bool bModify,void *ctx,
				uint32_t (*get_thread)  (void *,uint32_t,uint32_t,const void *),
				uint32_t (*hash_func)  (const void *key),
				int32_t (*get_tile)(void *ctx,const void *data,uint32_t nkey, void * const *keys,void const **grid),
				uint32_t pack_size,
				pack_function_type   pack_function,
				unpack_function_type unpack_function,
				init_function_type   init_function,
				uint32_t             flush_size,
				pack_function_type   flush_function,
				unpack_function_type combine_function,
				void* (*create_function) (void *,uint32_t,const void *))
	: CACHEhelper(nData,bModify), nPack(pack_size), nFlush(flush_size), ctx(ctx),
		get_thread(get_thread),hash_func(hash_func),get_tile(get_tile),
		pack_function(pack_function),unpack_function(unpack_function),
		init_function(init_function),flush_function(flush_function),
		combine_function(combine_function),create_function(create_function) {}
    };

extern "C"
void ramses_cache_open(MDL mdl,int cid,void *pHash,uint32_t (*hash_func)(const void *key),
				int iDataSize,bool bModify,void *ctx,
				uint32_t (*get_thread)  (void *,uint32_t,uint32_t,const void *),
				int32_t (*get_tile)(void *ctx,const void *data,uint32_t nkey, void * const *keys,void const **grid),
				uint32_t pack_size,
				ramsesCACHEhelper::  pack_function_type   pack_function,
				ramsesCACHEhelper::unpack_function_type unpack_function,
				ramsesCACHEhelper::  init_function_type   init_function,
				uint32_t flush_size,
				ramsesCACHEhelper::  pack_function_type   flush_function,
				ramsesCACHEhelper::unpack_function_type combine_function,
				void* (*create_function) (void *,uint32_t,const void *)) {
    auto hash = reinterpret_cast<hash::GHASH*>(pHash);
    reinterpret_cast<mdl::mdlClass *>(mdl)->AdvancedCacheInitialize(cid,hash,iDataSize,
		std::make_shared<ramsesCACHEhelper>(iDataSize,bModify,ctx,get_thread,hash_func,get_tile,
		pack_size,pack_function,unpack_function,init_function,
		flush_size,flush_function,combine_function,create_function));
    }
