# generated by datamodel-codegen:
#   filename:  https://opensource.ieee.org/2791-object/ieee-2791-schema/-/raw/master/io_domain.json
#   timestamp: 2022-09-13T23:51:51+00:00

from __future__ import annotations

from datetime import datetime
from enum import Enum
from typing import (
    List,
    Optional,
)

from pydantic import (
    AnyUrl,
    BaseModel,
    ConfigDict,
    EmailStr,
    Field,
    RootModel,
)


class Uri(BaseModel):
    model_config = ConfigDict(extra="forbid")

    filename: Optional[str] = None
    uri: AnyUrl
    access_time: Optional[datetime] = Field(
        None, description="Time stamp of when the request for this data was submitted"
    )
    sha1_checksum: Optional[str] = Field(
        None, description="output of hash function that produces a message digest", pattern="[A-Za-z0-9]+"
    )


class ObjectId(RootModel):
    root: str = Field(
        ...,
        description="A unique identifier that should be applied to each IEEE-2791 Object instance, generated and assigned by a IEEE-2791 database engine. IDs should never be reused",
    )


class ContributionEnum(Enum):
    authoredBy = "authoredBy"
    contributedBy = "contributedBy"
    createdAt = "createdAt"
    createdBy = "createdBy"
    createdWith = "createdWith"
    curatedBy = "curatedBy"
    derivedFrom = "derivedFrom"
    importedBy = "importedBy"
    importedFrom = "importedFrom"
    providedBy = "providedBy"
    retrievedBy = "retrievedBy"
    retrievedFrom = "retrievedFrom"
    sourceAccessedBy = "sourceAccessedBy"


class Contributor(BaseModel):
    model_config = ConfigDict(extra="forbid")

    name: str = Field(..., description="Name of contributor", examples=["Charles Darwin"])
    affiliation: Optional[str] = Field(
        None, description="Organization the particular contributor is affiliated with", examples=["HMS Beagle"]
    )
    email: Optional[EmailStr] = Field(
        None,
        description="electronic means for identification and communication purposes",
        examples=["name@example.edu"],
    )
    contribution: List[ContributionEnum] = Field(
        ..., description="type of contribution determined according to PAV ontology"
    )
    orcid: Optional[AnyUrl] = Field(
        None,
        description="Field to record author information. ORCID identifiers allow for the author to curate their information after submission. ORCID identifiers must be valid and must have the prefix ‘https://orcid.org/’",
        examples=["http://orcid.org/0000-0002-1825-0097"],
    )


class InputSubdomainItem(BaseModel):
    model_config = ConfigDict(extra="forbid")

    uri: Uri


class OutputSubdomainItem(BaseModel):
    mediatype: str = Field(
        ...,
        description="https://www.iana.org/assignments/media-types/",
        examples=["text/csv"],
        pattern="^(.*)$",
        title="mediatype",
    )
    uri: Uri


class InputAndOutputDomain(BaseModel):
    input_subdomain: List[InputSubdomainItem] = Field(
        ...,
        description="A record of the references and input files for the entire pipeline. Each type of input file is listed under a key for that type.",
        title="input_domain",
    )
    output_subdomain: List[OutputSubdomainItem] = Field(
        ..., description="A record of the outputs for the entire pipeline.", title="output_subdomain"
    )
