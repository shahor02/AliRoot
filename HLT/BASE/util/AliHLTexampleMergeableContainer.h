#ifndef AliHLTexampleMergeableContainer_h
#define AliHLTexampleMergeableContainer_h

#include "AliMergeable.h"
#include "AliHLTObjArray.h"
#include "TObject.h"

class AliHLTexampleMergeableContainer : public TObject, public AliMergeable {

private:
  AliHLTObjArray* fContainer;

public:
  AliHLTexampleMergeableContainer(): TObject(), AliMergeable(), fContainer(NULL) {
    //printf("AliHLTexampleMergeableContainer ctor, %p\n",this);
  }
  AliHLTexampleMergeableContainer(const char* name): TObject(), AliMergeable(),
    fContainer(new AliHLTObjArray(1))
  {
    fContainer->SetName(name);
    fContainer->SetOwner(kTRUE);
  }
  virtual ~AliHLTexampleMergeableContainer() {
    //printf("AliHLTexampleMergeableContainer dtor, %p\n",this);
    if (fContainer) fContainer->Delete();
    delete fContainer;
  }

  void Add(TObject*);

  virtual Long64_t Merge(TCollection *list);
  virtual TCollection* GetListOfDrawableObjects();

  ClassDef(AliHLTexampleMergeableContainer,1)
};
#endif
